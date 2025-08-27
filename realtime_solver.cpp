#include "hu.cpp" // Includes Card, HUState, HUGame, PokerEvaluator, action_to_string etc.
#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <sstream>
#include <random>
#include <algorithm>
#include <map>
#include <iomanip> // For std::fixed, std::setprecision if needed for amounts
#include <limits>  // Required for std::numeric_limits
#include <fstream>      // For preloadClusters
#include <mutex>        // For G_CLUSTER_CACHE_MUTEX, and now CFR mutexes
#include <unordered_map> // For G_CLUSTER_CACHE
#include <cctype>       // For ::isspace
#include <thread>       // For std::thread
#include <atomic>       // For atomic progress counter
#include <chrono>       // For timing in accelerated MCCFR

// Forward declarations for exploitability functions
double compute_current_strategy_exploitability(const HUState& start_state, int num_simulations = 1000);
double compute_strategy_file_exploitability(const HUState& start_state, 
                                           const std::string& strategy_csv_filename,
                                           int num_simulations = 1000);

// Forward declaration for accelerated MCCFR function
void solve_from_state_mccfr_accelerated(const HUState& start_state, 
                                       int num_iterations, 
                                       const std::string& output_csv_filename, 
                                       int num_threads_to_use,
                                       bool use_outcome_sampling,
                                       bool use_linear_cfr,
                                       int exploitability_check_interval);

// Global variables for clustering (mirrored from train.cpp)
std::unordered_map<std::string, std::unordered_map<std::string, std::string>> G_CLUSTER_CACHE;
std::mutex G_CLUSTER_CACHE_MUTEX;
bool G_CLUSTERS_LOADED = false;

// Global variables for CFR - optimized data structures
// These are cleared/reused by solve_from_state_cfr.
// Using unordered_map for O(1) average lookup instead of O(log n)
std::unordered_map<std::string, std::vector<double>> G_CFR_REGRET_SUM;
std::unordered_map<std::string, std::vector<double>> G_CFR_STRATEGY_SUM;
std::unordered_map<std::string, std::vector<std::pair<Action, double>>> G_CFR_INFOSET_ACTIONS; // Stores legal actions for an infoset

// Mutexes for CFR data structures
std::mutex G_CFR_REGRET_SUM_MUTEX;
std::mutex G_CFR_STRATEGY_SUM_MUTEX;
std::mutex G_CFR_INFOSET_ACTIONS_MUTEX; // Guards G_CFR_INFOSET_ACTIONS and initial population of regret/strategy sums for new infosets

std::atomic<int> G_CFR_COMPLETED_ITERATIONS_COUNT(0); // For progress reporting

// New: Linear CFR iteration tracking
std::atomic<int> G_CFR_CURRENT_ITERATION(0); // For Linear CFR weighting
bool G_USE_LINEAR_CFR = true; // Enable Linear CFR by default for faster convergence

// New: Subgame CFR globals
std::mt19937 G_SUBGAME_RNG(std::random_device{}()); // RNG for subgame sampling
std::mutex G_SUBGAME_RNG_MUTEX; // Mutex for RNG access
std::vector<Card> G_ORIGINAL_COMMUNITY_CARDS; // Store the original community cards input by user

// If PokerEvaluator in hu.cpp doesn't have RANK_ORDER_FOR_PATTERN, this can be used.
// train.cpp uses evaluator.RANK_ORDER_FOR_PATTERN, implying it's a member.
// If your PokerEvaluator has it, you can adapt get_rank_pattern_string_from_sorted_hand_ic
// to use evaluator.RANK_ORDER_FOR_PATTERN.
const std::vector<std::string> RANK_ORDER_FOR_PATTERN_CONST_IC = {
    "2", "3", "4", "5", "6", "7", "8", "9", "10", "J", "Q", "K", "A"
};

// Helper function to sort a poker hand by rank and return it as a string
std::string sort_cards_by_rank_to_string_ic(const std::vector<Card>& cards) {
    std::vector<Card> sortedCards = cards;
    PokerEvaluator evaluator; // For RANK_VALUES

    std::sort(sortedCards.begin(), sortedCards.end(), [&](const Card& a, const Card& b) {
        return evaluator.RANK_VALUES.at(a.rank) < evaluator.RANK_VALUES.at(b.rank);
    });

    std::string sortedHandStr;
    for (size_t i = 0; i < sortedCards.size(); ++i) {
        sortedHandStr += sortedCards[i].toString();
        if (i < sortedCards.size() - 1) {
            sortedHandStr += " ";
        }
    }
    return sortedHandStr;
}

// Helper function to parse a sorted hand string and return a rank pattern string
// e.g., "Ac Kd Qh" -> "{'A':1,'K':1,'Q':1}" (simplified, actual output is more specific)
std::string get_rank_pattern_string_from_sorted_hand_ic(const std::string& sorted_hand_str) {
    std::vector<std::string> cards_str_vec;
    std::istringstream iss(sorted_hand_str);
    std::string card_str;
    
    while (iss >> card_str) {
        cards_str_vec.push_back(card_str);
    }

    std::map<std::string, int> rank_counts;
    for (const auto& c_str : cards_str_vec) {
        if (!c_str.empty()) {
            std::string rank_val_str;
            if (c_str.length() > 1) { // Handles "10s" -> "10", "As" -> "A"
                rank_val_str = c_str.substr(0, c_str.length() - 1);
            } else { 
                rank_val_str = c_str; // Fallback, though unlikely for valid card strings
            }
            rank_counts[rank_val_str]++;
        }
    }

    std::ostringstream result;
    result << "\"{\'"; // e.g., "{'
    bool first = true;
    for (const auto& rank_key : RANK_ORDER_FOR_PATTERN_CONST_IC) {
        if (rank_counts.count(rank_key) > 0) {
            if (!first) {
                result << ", \'"; // e.g., , '
            }
            result << rank_key << "\':" << rank_counts[rank_key]; // MODIFICATION: No space after colon. e.g., K':1
            first = false;
        }
    }
    result << "}\""; // e.g., }"

    return result.str(); // Expected: "{'A':1,'K':1}"
}

// Helper function to get suit pattern string from a sorted hand string
// This is a complex function adapted from inspiration.cpp AND train2.cpp logic
std::string get_suit_pattern_string_from_sorted_hand_ic(const std::string& sorted_hand_str) {
    std::vector<std::string> card_strings;
    std::istringstream iss(sorted_hand_str);
    std::string card_token;
    while (iss >> card_token) {
        card_strings.push_back(card_token);
    }
    
    std::vector<std::string> hand_ranks;
    std::vector<char> hand_suits;
    PokerEvaluator evaluator; // For RANK_VALUES

    for (const auto& c_str : card_strings) {
        if (c_str.length() > 1) {
            hand_ranks.push_back(c_str.substr(0, c_str.length() - 1));
            hand_suits.push_back(c_str.back());
        } else if (!c_str.empty()) {
            hand_ranks.push_back(c_str);
            // Assuming valid card strings, suit should always be present.
            // If suit can be missing, error handling or default suit needed.
        }
    }
    
    // Create rank-to-group mapping (inspired by train2.cpp's suit_pattern)
    std::unordered_map<std::string, int> rank_to_group_mapping;
    int group_idx_assigner = 0;
    
    std::map<std::string, int> current_hand_rank_counts;
    for (const auto& r : hand_ranks) {
        current_hand_rank_counts[r]++;
    }
    
    std::vector<std::pair<std::string, int>> sorted_ranks_with_counts;
    for (const auto& pair_entry : current_hand_rank_counts) {
        sorted_ranks_with_counts.push_back(pair_entry);
    }
    
    // Sort by count (descending), then by rank value (descending using PokerEvaluator)
    std::sort(sorted_ranks_with_counts.begin(), sorted_ranks_with_counts.end(),
        [&evaluator](const auto& a, const auto& b) {
            if (a.second != b.second) {
                return a.second > b.second; // Sort by count (descending)
            }
            // Fallback to rank value for tie-breaking (descending)
            return evaluator.RANK_VALUES.at(a.first) > evaluator.RANK_VALUES.at(b.first);
        });
    
    for (const auto& pair_entry : sorted_ranks_with_counts) {
        rank_to_group_mapping[pair_entry.first] = group_idx_assigner++;
    }
    
    // Create signature for each suit
    // Using std::map for suit_signature to ensure suits are processed in a consistent order (c,d,h,s)
    // which helps in the final sorting of the pattern if suits themselves are part of the sort criteria.
    // train2.cpp uses unordered_map then sorts the resulting vector of vectors.
    std::map<char, std::vector<int>> suit_signatures;
    for (size_t i = 0; i < card_strings.size(); ++i) {
        if (i >= hand_suits.size() || i >= hand_ranks.size()) continue; // Safety check

        char current_suit = hand_suits[i];
        std::string current_rank_str = hand_ranks[i];
        
        // Ensure rank exists in mapping (it should if logic is correct)
        if (rank_to_group_mapping.find(current_rank_str) == rank_to_group_mapping.end()) {
            // This would indicate an issue, perhaps a rank in hand_ranks not in current_hand_rank_counts
            // or an empty hand_ranks/hand_suits for a non-empty card_strings.
            // For robustness, one might add error logging here.
            continue;
        }
        int group_val = rank_to_group_mapping.at(current_rank_str);
        
        // Add group_val to suit's signature only if not already present
        // This handles multiple cards of the same "group" for a given suit (e.g. from community cards)
        // correctly reflecting the distinct groups present for that suit.
        bool already_in_signature = false;
        for(int sig_val : suit_signatures[current_suit]) {
            if (sig_val == group_val) {
                already_in_signature = true;
                break;
            }
        }
        if (!already_in_signature) {
            suit_signatures[current_suit].push_back(group_val);
        }
    }
    
    // Sort each suit's signature (the group indices within each suit's vector)
    for (auto& sig_pair : suit_signatures) {
        std::sort(sig_pair.second.begin(), sig_pair.second.end());
    }
    
    // Convert to vector of vectors (pattern) and sort for consistent output
    // The outer sort is by the suit signatures themselves (e.g. [0,1] comes before [0,2])
    std::vector<std::vector<int>> final_pattern;
    for (const auto& sig_pair : suit_signatures) {
        if (!sig_pair.second.empty()) { // Only add if the suit had cards contributing to a group
            final_pattern.push_back(sig_pair.second);
        }
    }
    std::sort(final_pattern.begin(), final_pattern.end());
    
    // Format the output string
    std::ostringstream result_oss;
    result_oss << "\"["; // e.g., "[
    
    for (size_t i = 0; i < final_pattern.size(); ++i) {
        result_oss << "[";
        for (size_t j = 0; j < final_pattern[i].size(); ++j) {
            result_oss << final_pattern[i][j];
            if (j < final_pattern[i].size() - 1) {
                result_oss << ","; // MODIFICATION: No space after comma
            }
        }
        result_oss << "]";
        if (i < final_pattern.size() - 1) {
            result_oss << ","; // MODIFICATION: No space after comma
        }
    }
    
    result_oss << "]\""; // e.g., ]"
    return result_oss.str(); // Expected: "[[0,1],[2]]"
}

// Function to preload all clusters from CSV files (from train.cpp)
void preloadClusters_ic() {
    std::lock_guard<std::mutex> lock(G_CLUSTER_CACHE_MUTEX);
    if (G_CLUSTERS_LOADED) return;

    std::cout << "Preloading clusters for infoset_counter..." << std::endl;
    std::vector<std::string> rounds = {"flop", "turn", "river"};
    // Path to cluster files, ensure these exist and are correctly formatted.
    // These paths are relative to the execution directory of infoset_counter.
    std::string base_path = "./utils/"; // Or make this configurable

    for (const auto& round : rounds) {
        std::string filename = base_path + "repr_l2_" + round + "_equities_clustered_hands_with_avg_equity.csv";
        std::ifstream file(filename);

        if (!file.is_open()) {
            std::cerr << "Error: Could not open cluster file " << filename << std::endl;
            continue;
        }

        G_CLUSTER_CACHE[round].reserve(10000); 
        std::string line;
        std::getline(file, line); // Skip header

        while (std::getline(file, line)) {
            std::vector<std::string> fields;
            std::string field;
            bool inQuotes = false;
            // Simple CSV parser from train.cpp
            for (char c : line) {
                if (c == '"') {
                    inQuotes = !inQuotes;
                    field += c; // Keep quotes as part of the field
                } else if (c == ',' && !inQuotes) {
                    fields.push_back(field);
                    field.clear();
                } else {
                    field += c;
                }
            }
            fields.push_back(field); // Add the last field


            if (fields.size() >= 8) { // Ensure enough fields (cluster_id=fields[3], rank_pattern=fields[6], suit_pattern=fields[7])
                std::string cluster_str = fields[3];
                std::string rank_pattern_csv = fields[6];
                std::string suit_pattern_csv = fields[7];
                
                // Clean whitespace from patterns read from CSV, similar to train2.cpp
                std::string rank_pattern_key = rank_pattern_csv;
                rank_pattern_key.erase(std::remove_if(rank_pattern_key.begin(), rank_pattern_key.end(), ::isspace), rank_pattern_key.end());
                
                std::string suit_pattern_key = suit_pattern_csv;
                suit_pattern_key.erase(std::remove_if(suit_pattern_key.begin(), suit_pattern_key.end(), ::isspace), suit_pattern_key.end());

                std::string key = rank_pattern_key + "|" + suit_pattern_key;
                G_CLUSTER_CACHE[round][key] = cluster_str;
            }
        }
        std::cout << "Loaded " << G_CLUSTER_CACHE[round].size() << " clusters for " << round << " from " << filename << std::endl;
        file.close();
    }
    G_CLUSTERS_LOADED = true;
    std::cout << "Cluster loading complete for infoset_counter." << std::endl;
}

// Updated function to find cluster using the cache
std::string find_cluster_from_cache_ic(const std::string& round, const std::string& rankPattern, const std::string& suitPattern) {
    if (!G_CLUSTERS_LOADED) {
        std::cerr << "Warning: Clusters not preloaded. Attempting to load now." << std::endl;
        preloadClusters_ic();
        if (!G_CLUSTERS_LOADED) {
            std::cerr << "Error: Cluster loading failed." << std::endl;
            return "CLUST_LOAD_ERR";
        }
    }

    if (round == "preflop") {
        return "N/A";
    }

    std::string cleanRankPattern = rankPattern;
    cleanRankPattern.erase(std::remove_if(cleanRankPattern.begin(), cleanRankPattern.end(), ::isspace), cleanRankPattern.end());

    std::string cleanSuitPattern = suitPattern;
    cleanSuitPattern.erase(std::remove_if(cleanSuitPattern.begin(), cleanSuitPattern.end(), ::isspace), cleanSuitPattern.end());

    std::string key = cleanRankPattern + "|" + cleanSuitPattern;
    
    std::lock_guard<std::mutex> lock(G_CLUSTER_CACHE_MUTEX);
    
    if (!G_CLUSTER_CACHE.count(round)) {
        std::cerr << "Round '" << round << "' not found in cluster cache!" << std::endl;
        return "NOT_FOUND_ROUND_MISSING";
    }

    if (G_CLUSTER_CACHE.at(round).count(key)) {
        return G_CLUSTER_CACHE.at(round).at(key);
    } else {
        std::cerr << "Key not found in " << round << " cache!" << std::endl;
        return "NOT_FOUND_KEY_MISSING";
    }
}

// Get board cluster ID, adapted for infoset_counter (uses int round from HUState)
std::string getBoardClusterID_ic(int game_round_int, // 0=preflop, 1=flop, 2=turn, 3=river
                              const std::vector<Card>& community_cards,
                              const Card& h1, const Card& h2) {
    std::string round_name;
    std::vector<Card> relevant_cards;
    relevant_cards.push_back(h1);
    relevant_cards.push_back(h2);

    int num_community_for_round = 0;

    switch (game_round_int) {
        case 1: // Flop
            round_name = "flop";
            num_community_for_round = 3;
            break;
        case 2: // Turn
            round_name = "turn";
            num_community_for_round = 4;
            break;
        case 3: // River
            round_name = "river";
            num_community_for_round = 5;
            break;
        default:
            return "N/A_ROUND"; // Not a post-flop round for this clustering or invalid round
    }

    if (community_cards.size() < static_cast<size_t>(num_community_for_round)) {
        // std::cerr << "Warning: Not enough community cards for " << round_name << " cluster. Have " << community_cards.size() << ", need " << num_community_for_round << std::endl;
        return "N/A_CARDS";
    }

    // Add all relevant community cards
    for (int i = 0; i < num_community_for_round; ++i) {
        relevant_cards.push_back(community_cards[i]);
    }

    // Debug output
    //std::cout << "DEBUG: Processing " << round_name << " with " << relevant_cards.size() << " cards" << std::endl;
    
    std::string sorted_hand_str = sort_cards_by_rank_to_string_ic(relevant_cards);
    //std::cout << "DEBUG: Sorted hand string: " << sorted_hand_str << std::endl;
    
    std::string rank_pattern = get_rank_pattern_string_from_sorted_hand_ic(sorted_hand_str);
    std::string suit_pattern = get_suit_pattern_string_from_sorted_hand_ic(sorted_hand_str);
    
    //std::cout << "DEBUG: Generated patterns for lookup:" << std::endl;
    //std::cout << "  Rank pattern: " << rank_pattern << std::endl;
    //std::cout << "  Suit pattern: " << suit_pattern << std::endl;
    
    return find_cluster_from_cache_ic(round_name, rank_pattern, suit_pattern);
}

// Helper function to get preflop hand abstraction (e.g., AKs, 77, T9o)
std::string get_preflop_abstraction_ic(const Card& c1, const Card& c2) {
    PokerEvaluator evaluator; // For RANK_VALUES

    std::string r1_val = c1.rank;
    std::string s1_val = c1.suit;
    std::string r2_val = c2.rank;
    std::string s2_val = c2.suit;

    // Ensure r1_val is the higher or equal rank for consistent ordering (e.g., AKs not KAs)
    // The PokerEvaluator::RANK_VALUES map should map ranks to comparable integers.
    if (evaluator.RANK_VALUES.at(r1_val) < evaluator.RANK_VALUES.at(r2_val)) {
        std::swap(r1_val, r2_val);
        std::swap(s1_val, s2_val); // Keep suit paired with its rank
    }

    std::string abstraction;
    if (r1_val == r2_val) { // Pair
        abstraction = r1_val + r2_val;
    } else { // Non-pair
        abstraction = r1_val + r2_val;
        if (s1_val == s2_val) {
            abstraction += "s"; // Suited
        } else {
            abstraction += "o"; // Offsuit
        }
    }
    return abstraction;
}

// Helper function to convert round string to int
int getRoundInt(const std::string& round_str) {
    if (round_str == "preflop") return 0;
    if (round_str == "flop") return 1;
    if (round_str == "turn") return 2;
    if (round_str == "river") return 3;
    return -1; // Invalid round
}

// Generates a string key for the current player's information set
std::string get_infoset_key_for_state(const HUState& state) {
    std::ostringstream oss;

    // Round
    oss << state.round << "|";

    // Current player
    oss << state.current_player() << "|";

    // Last bet amount
    oss << state.current_bet_to_match << "|";

    // Cumulated pot
    oss << state.cumulative_pot << "|";

    int round_int = getRoundInt(state.round);
    int player_idx = state.current_player();

    // Add either preflop abstraction or postflop cluster
    if (round_int == 0) { // Preflop
        if (state.cards.size() >= (player_idx * 2) + 2) {
            Card p_h1 = state.cards[player_idx * 2];
            Card p_h2 = state.cards[player_idx * 2 + 1];
            oss << get_preflop_abstraction_ic(p_h1, p_h2);
        } else {
            oss << "N/A_CARDS";
        }
    } else if (round_int >= 1 && round_int <= 3) { // Postflop (Flop, Turn, River)
        if (state.cards.size() >= (player_idx * 2) + 2 && !state.community_cards.empty()) {
            Card p_h1 = state.cards[player_idx * 2];
            Card p_h2 = state.cards[player_idx * 2 + 1];
            oss << getBoardClusterID_ic(round_int, state.community_cards, p_h1, p_h2);
        } else {
            oss << "N/A_INFOS";
        }
    }

    return oss.str();
}

// Simulates games using Monte Carlo and counts unique information sets
// MODIFIED: Now starts from a given state_to_simulate_from
int count_unique_infosets_mc(
    const HUState& state_to_simulate_from, // The state to start simulations from
    int num_simulations_from_node,
    int max_actions_per_simulation) {

    std::set<std::string> unique_infosets;
    std::mt19937 rng(std::random_device{}()); // Random number generator

    std::cout << "Starting Monte Carlo simulation for " << num_simulations_from_node
              << " games from the current state..." << std::endl;

    for (int i = 0; i < num_simulations_from_node; ++i) {
        if ((i + 1) % (num_simulations_from_node / 10 == 0 ? 1 : num_simulations_from_node / 10) == 0) {
            std::cout << "Sub-simulation " << (i + 1) << "/" << num_simulations_from_node
                      << ", Unique Infosets found so far: " << unique_infosets.size() << std::endl;
        }

        HUState current_sim_state = state_to_simulate_from; // Make a copy to simulate from

        for (int action_count = 0; action_count < max_actions_per_simulation; ++action_count) {
            if (current_sim_state.game_over) {
                break;
            }

            if (current_sim_state.is_chance_node()) {
                // Optional: could generate an infoset key for chance nodes if desired
                // std::string chance_key = get_infoset_key_for_state(current_sim_state);
                // unique_infosets.insert(chance_key);
                current_sim_state.apply_action(Action::DEAL);
            } else {
                // Player node
                std::string infoset_key = get_infoset_key_for_state(current_sim_state);
                unique_infosets.insert(infoset_key);

                std::vector<std::pair<Action, double>> legal_actions = current_sim_state.legal_actions(
                    current_sim_state.game.pot_fraction_bet_sizes,
                    current_sim_state.game.fixed_bet_sizes_bb
                );
                if (legal_actions.empty()) {
                    // This should ideally not happen if game is not over and not a chance node.
                    break;
                }

                std::uniform_int_distribution<size_t> dist(0, legal_actions.size() - 1);
                size_t chosen_action_idx = dist(rng);
                std::pair<Action, double> chosen_action = legal_actions[chosen_action_idx];

                current_sim_state.apply_action(chosen_action.first, chosen_action.second);
            }
        }
         // Optional: Add infoset for terminal state if game ended
        if(current_sim_state.game_over){
            // std::string terminal_key = get_infoset_key_for_state(current_sim_state);
            // unique_infosets.insert(terminal_key);
        }
    }
    std::cout << "Monte Carlo sub-simulation finished." << std::endl;
    return unique_infosets.size();
}

// MCCFR Helper: Get current strategy based on regrets (Regret Matching) - CFR+ Version
std::vector<double> get_mccfr_strategy(const std::string& infoset_key, int num_actions) {
    std::vector<double> strategy(num_actions);
    double positive_regret_sum = 0.0;

    { // Scope for the lock
        std::lock_guard<std::mutex> lock(G_CFR_REGRET_SUM_MUTEX);
        // Assumes infoset_key exists in G_CFR_REGRET_SUM due to prior initialization.
        const auto& regrets = G_CFR_REGRET_SUM.at(infoset_key); 
        for (int a = 0; a < num_actions; ++a) {
            // CFR+: Use max(regret, 0) for strategy computation
            strategy[a] = std::max(regrets[a], 0.0);
            positive_regret_sum += strategy[a];
        }
    } // G_CFR_REGRET_SUM_MUTEX is released here

    if (positive_regret_sum > 0) {
        for (int a = 0; a < num_actions; ++a) {
            strategy[a] /= positive_regret_sum;
        }
    } else {
        // Default to uniform random if all regrets are non-positive
        for (int a = 0; a < num_actions; ++a) {
            strategy[a] = 1.0 / static_cast<double>(num_actions);
        }
    }
    return strategy;
}

// MCCFR Recursive Function using External Sampling
// Returns expected utility for the player being updated
double mccfr_recursive(HUState current_state, int updating_player, std::mt19937& rng) {
    if (current_state.game_over) {
        std::vector<double> returns = current_state.returns();
        return returns[updating_player];
    }

    if (current_state.is_chance_node()) {
        HUState next_state = current_state;
        next_state.apply_action(Action::DEAL);
        return mccfr_recursive(next_state, updating_player, rng);
    }

    int current_player = current_state.current_player();
    std::string infoset_key = get_infoset_key_for_state(current_state);

    std::vector<std::pair<Action, double>> legal_actions;
    int num_actions_for_infoset;

    { // Scope for G_CFR_INFOSET_ACTIONS_MUTEX
        std::lock_guard<std::mutex> lock(G_CFR_INFOSET_ACTIONS_MUTEX);
        auto it_actions = G_CFR_INFOSET_ACTIONS.find(infoset_key);
        if (it_actions == G_CFR_INFOSET_ACTIONS.end()) {
            legal_actions = current_state.legal_actions(
                current_state.game.pot_fraction_bet_sizes,
                current_state.game.fixed_bet_sizes_bb
            );
            if (legal_actions.empty()) {
                std::cerr << "MCCFR Warning: Player " << current_player << " has no legal actions at infoset "
                          << infoset_key << " but game not over and not chance node. State:\n"
                          << current_state.to_string(false) << std::endl;
                std::vector<double> returns = current_state.returns();
                return returns[updating_player];
            }
            G_CFR_INFOSET_ACTIONS[infoset_key] = legal_actions;
            G_CFR_REGRET_SUM.try_emplace(infoset_key, std::vector<double>(legal_actions.size(), 0.0));
            G_CFR_STRATEGY_SUM.try_emplace(infoset_key, std::vector<double>(legal_actions.size(), 0.0));
            num_actions_for_infoset = legal_actions.size();
        } else {
            legal_actions = it_actions->second;
            num_actions_for_infoset = legal_actions.size();
        }
    }

    if (num_actions_for_infoset == 0) { 
        std::cerr << "MCCFR Error: num_actions is 0 for infoset " << infoset_key << std::endl;
        std::vector<double> returns = current_state.returns();
        return returns[updating_player];
    }

    std::vector<double> strategy = get_mccfr_strategy(infoset_key, num_actions_for_infoset);

    if (current_player == updating_player) {
        // Player being updated: traverse all actions
        std::vector<double> action_utilities(num_actions_for_infoset);
        double node_utility = 0.0;

        for (int a = 0; a < num_actions_for_infoset; ++a) {
            HUState next_state = current_state;
            next_state.apply_action(legal_actions[a].first, legal_actions[a].second);
            
            action_utilities[a] = mccfr_recursive(next_state, updating_player, rng);
            node_utility += strategy[a] * action_utilities[a];
        }

        // Update regrets
        { 
            std::lock_guard<std::mutex> lock(G_CFR_REGRET_SUM_MUTEX);
            
            // Ensure infoset exists in regret sum
            if (G_CFR_REGRET_SUM.find(infoset_key) == G_CFR_REGRET_SUM.end()) {
                G_CFR_REGRET_SUM[infoset_key] = std::vector<double>(num_actions_for_infoset, 0.0);
            }
            
            // Ensure correct size
            if (static_cast<int>(G_CFR_REGRET_SUM[infoset_key].size()) != num_actions_for_infoset) {
                G_CFR_REGRET_SUM[infoset_key].resize(num_actions_for_infoset, 0.0);
            }
            
            for (int a = 0; a < num_actions_for_infoset; ++a) {
                double regret = action_utilities[a] - node_utility;
                // CFR+: Clamp regrets to be non-negative
                G_CFR_REGRET_SUM[infoset_key][a] = std::max(G_CFR_REGRET_SUM[infoset_key][a] + regret, 0.0);
            }
        }
        
        // Update average strategy
        { 
            std::lock_guard<std::mutex> lock(G_CFR_STRATEGY_SUM_MUTEX);
            
            // Ensure infoset exists in strategy sum
            if (G_CFR_STRATEGY_SUM.find(infoset_key) == G_CFR_STRATEGY_SUM.end()) {
                G_CFR_STRATEGY_SUM[infoset_key] = std::vector<double>(num_actions_for_infoset, 0.0);
            }
            
            // Ensure correct size
            if (static_cast<int>(G_CFR_STRATEGY_SUM[infoset_key].size()) != num_actions_for_infoset) {
                G_CFR_STRATEGY_SUM[infoset_key].resize(num_actions_for_infoset, 0.0);
            }
            
            // Linear CFR: Weight strategy updates by iteration number for faster convergence
            double weight = G_USE_LINEAR_CFR ? std::max(1, G_CFR_CURRENT_ITERATION.load()) : 1.0;
            
            for (int a = 0; a < num_actions_for_infoset; ++a) {
                G_CFR_STRATEGY_SUM[infoset_key][a] += weight * strategy[a];
            }
        }
        
        return node_utility;
    } else {
        // Opponent: sample action according to strategy
        std::discrete_distribution<int> action_dist(strategy.begin(), strategy.end());
        int sampled_action = action_dist(rng);
        
        HUState next_state = current_state;
        next_state.apply_action(legal_actions[sampled_action].first, legal_actions[sampled_action].second);
        
        return mccfr_recursive(next_state, updating_player, rng);
    }
}

// Main function to run MCCFR from a given state
void solve_from_state_mccfr(const HUState& start_state, int num_iterations, const std::string& output_csv_filename, int num_threads_to_use) {
    // Clear global CFR data structures for a new solve
    G_CFR_REGRET_SUM.clear();
    G_CFR_STRATEGY_SUM.clear();
    G_CFR_INFOSET_ACTIONS.clear();
    G_CFR_COMPLETED_ITERATIONS_COUNT = 0;

    std::cout << "Starting MCCFR for " << num_iterations << " iterations using " << num_threads_to_use << " threads." << std::endl;
    std::cout << "Output will be saved to: " << output_csv_filename << std::endl;

    std::vector<std::thread> threads;
    int iterations_per_thread_base = (num_iterations > 0 && num_threads_to_use > 0) ? (num_iterations / num_threads_to_use) : 0;
    int iterations_remainder = (num_iterations > 0 && num_threads_to_use > 0) ? (num_iterations % num_threads_to_use) : 0;

    auto mccfr_task = [&](int iterations_for_this_thread, int thread_id) {
        std::mt19937 thread_rng(std::random_device{}() + thread_id);
        
        for (int i = 0; i < iterations_for_this_thread; ++i) {
            // Increment iteration counter for Linear CFR weighting
            int current_iter = G_CFR_CURRENT_ITERATION.fetch_add(1, std::memory_order_relaxed);
            
            // Alternate between updating player 0 and player 1
            int updating_player = (G_CFR_COMPLETED_ITERATIONS_COUNT.load(std::memory_order_relaxed) + i) % 2;
            mccfr_recursive(start_state, updating_player, thread_rng);
            G_CFR_COMPLETED_ITERATIONS_COUNT.fetch_add(1, std::memory_order_relaxed);
        }
    };

    for (int t = 0; t < num_threads_to_use; ++t) {
        int iterations_for_this_thread = iterations_per_thread_base + (t < iterations_remainder ? 1 : 0);
        if (iterations_for_this_thread > 0) {
            threads.emplace_back(mccfr_task, iterations_for_this_thread, t);
        }
    }

    // Progress monitoring by main thread
    int report_interval_ms = 500;
    int iterations_reported_at_last_print = -1;
    bool first_report_triggered = false;

    if (num_iterations > 0) {
        std::cout << "MCCFR Iteration progress (target " << num_iterations << "):" << std::endl;
    }

    while(G_CFR_COMPLETED_ITERATIONS_COUNT.load(std::memory_order_acquire) < num_iterations) {
        std::this_thread::sleep_for(std::chrono::milliseconds(report_interval_ms));
        int current_completed = G_CFR_COMPLETED_ITERATIONS_COUNT.load(std::memory_order_acquire);
        
        if (!first_report_triggered || current_completed > iterations_reported_at_last_print) {
            std::cout << "  Completed " << current_completed << " / " << num_iterations << " iterations." << std::endl;
            iterations_reported_at_last_print = current_completed;
            first_report_triggered = true;
        }
        if (current_completed >= num_iterations) break; 
    }

    // Join threads
    for (auto& th : threads) {
        if (th.joinable()) {
            th.join();
        }
    }
    
    // Final completion message
    int final_completed_count = G_CFR_COMPLETED_ITERATIONS_COUNT.load(std::memory_order_acquire);
    if (num_iterations > 0) {
        if (final_completed_count > iterations_reported_at_last_print || !first_report_triggered) {
             std::cout << "  Completed " << final_completed_count << " / " << num_iterations << " iterations." << std::endl;
        }
        std::cout << "MCCFR Iteration " << final_completed_count << "/" << num_iterations << " (All threads complete)" << std::endl;
    }

    std::cout << "MCCFR iterations complete. Saving strategy..." << std::endl;
    std::ofstream outfile(output_csv_filename);
    outfile << "infoset_key,action_string,probability\n";

    for (const auto& pair : G_CFR_INFOSET_ACTIONS) {
        const std::string& infoset_key = pair.first;
        const auto& actions_for_infoset = pair.second;

        if (G_CFR_STRATEGY_SUM.count(infoset_key)) {
            const std::vector<double>& summed_strategy = G_CFR_STRATEGY_SUM.at(infoset_key);
            double total_summed_strategy = 0.0;
            for (double prob_sum : summed_strategy) {
                total_summed_strategy += prob_sum;
            }

            if (total_summed_strategy > 1e-9) {
                for (size_t a = 0; a < summed_strategy.size(); ++a) {
                    double avg_prob = summed_strategy[a] / total_summed_strategy;
                    outfile << "\"" << infoset_key << "\",\"" 
                            << action_to_string(actions_for_infoset[a].first, actions_for_infoset[a].second) 
                            << "\"," << std::fixed << std::setprecision(5) << avg_prob << "\n";
                }
            } else {
                if (!actions_for_infoset.empty()) {
                    double uniform_prob = 1.0 / static_cast<double>(actions_for_infoset.size());
                    for (size_t a = 0; a < actions_for_infoset.size(); ++a) {
                        outfile << "\"" << infoset_key << "\",\""
                                << action_to_string(actions_for_infoset[a].first, actions_for_infoset[a].second)
                                << "\"," << std::fixed << std::setprecision(5) << uniform_prob << "\n";
                    }
                }
            }
        } else {
            std::cerr << "Warning: Infoset " << infoset_key << " found in actions but not in strategy sum." << std::endl;
            if (!actions_for_infoset.empty()) {
                double uniform_prob = 1.0 / static_cast<double>(actions_for_infoset.size());
                 for (size_t a = 0; a < actions_for_infoset.size(); ++a) {
                    outfile << "\"" << infoset_key << "\",\""
                            << action_to_string(actions_for_infoset[a].first, actions_for_infoset[a].second)
                            << "\"," << std::fixed << std::setprecision(5) << uniform_prob << "\n";
                }
            }
        }
    }
    outfile.close();
    std::cout << "Strategy saved to " << output_csv_filename << std::endl;
}

// MCCFR solving with periodic exploitability computation for convergence monitoring
void solve_from_state_mccfr_with_exploitability_tracking(const HUState& start_state, 
                                                         int num_iterations, 
                                                         const std::string& output_csv_filename, 
                                                         int num_threads_to_use,
                                                         int exploitability_check_interval = 100,
                                                         int exploitability_simulations = 500) {
    // Clear global CFR data structures for a new solve
    G_CFR_REGRET_SUM.clear();
    G_CFR_STRATEGY_SUM.clear();
    G_CFR_INFOSET_ACTIONS.clear();
    G_CFR_COMPLETED_ITERATIONS_COUNT = 0;

    // Get starting pot for percentage calculation
    double starting_pot = start_state.cumulative_pot;
    std::cout << "Starting pot: " << starting_pot << " BB" << std::endl;

    std::cout << "Starting MCCFR with exploitability tracking for " << num_iterations 
              << " iterations using " << num_threads_to_use << " threads." << std::endl;
    std::cout << "Exploitability will be computed every " << exploitability_check_interval 
              << " iterations using " << exploitability_simulations << " simulations." << std::endl;
    std::cout << "Output will be saved to: " << output_csv_filename << std::endl;

    std::vector<std::thread> threads;
    int iterations_per_thread_base = (num_iterations > 0 && num_threads_to_use > 0) ? (num_iterations / num_threads_to_use) : 0;
    int iterations_remainder = (num_iterations > 0 && num_threads_to_use > 0) ? (num_iterations % num_threads_to_use) : 0;

    // Exploitability tracking variables
    std::vector<std::tuple<int, double, double>> exploitability_history; // (iteration, exploitability, exploitability_percentage)
    std::mutex exploitability_history_mutex;
    std::atomic<int> next_exploitability_check(exploitability_check_interval);

    auto mccfr_task = [&](int iterations_for_this_thread, int thread_id) {
        std::mt19937 thread_rng(std::random_device{}() + thread_id);
        
        for (int i = 0; i < iterations_for_this_thread; ++i) {
            // Increment iteration counter for Linear CFR weighting
            int current_iter = G_CFR_CURRENT_ITERATION.fetch_add(1, std::memory_order_relaxed);
            
            // Alternate between updating player 0 and player 1
            int updating_player = (G_CFR_COMPLETED_ITERATIONS_COUNT.load(std::memory_order_relaxed) + i) % 2;
            mccfr_recursive(start_state, updating_player, thread_rng);
            
            int completed_count = G_CFR_COMPLETED_ITERATIONS_COUNT.fetch_add(1, std::memory_order_relaxed) + 1;
            
            // Check if this thread should perform an exploitability check
            int expected_check = next_exploitability_check.load(std::memory_order_acquire);
            if (completed_count >= expected_check && completed_count > 0) {
                // Try to claim this exploitability check
                if (next_exploitability_check.compare_exchange_strong(expected_check, expected_check + exploitability_check_interval)) {
                    std::cout << "\n--- Exploitability Check at Iteration " << completed_count << " (Thread " << thread_id << ") ---" << std::endl;
                    
                    try {
                        double exploitability = compute_current_strategy_exploitability(start_state, exploitability_simulations);
                        double exploitability_percentage = (starting_pot > 0) ? (exploitability / starting_pot) * 100.0 : 0.0;
                        
                        {
                            std::lock_guard<std::mutex> lock(exploitability_history_mutex);
                            exploitability_history.push_back(std::make_tuple(completed_count, exploitability, exploitability_percentage));
                        }
                        
                        std::cout << "Iteration " << completed_count << " - Exploitability: " << std::fixed 
                                  << std::setprecision(6) << exploitability << " BB/hand" << std::endl;
                        std::cout << "Iteration " << completed_count << " - Exploitability Percentage: " << std::fixed 
                                  << std::setprecision(2) << exploitability_percentage << "% of starting pot" << std::endl;
                        
                        // Show improvement if we have previous data
                        {
                            std::lock_guard<std::mutex> lock(exploitability_history_mutex);
                            if (exploitability_history.size() > 1) {
                                double prev_exploitability = std::get<1>(exploitability_history[exploitability_history.size() - 2]);
                                double prev_exploitability_percentage = std::get<2>(exploitability_history[exploitability_history.size() - 2]);
                                double improvement = prev_exploitability - exploitability;
                                double improvement_percentage = prev_exploitability_percentage - exploitability_percentage;
                                std::cout << "Improvement since last check: " << std::fixed 
                                          << std::setprecision(6) << improvement << " BB/hand (" 
                                          << std::setprecision(2) << improvement_percentage << " percentage points)" << std::endl;
                            }
                        }
                        std::cout << "--- End Exploitability Check ---\n" << std::endl;
                        
                    } catch (const std::exception& e) {
                        std::cerr << "Error computing exploitability: " << e.what() << std::endl;
                    }
                }
            }
        }
    };

    for (int t = 0; t < num_threads_to_use; ++t) {
        int iterations_for_this_thread = iterations_per_thread_base + (t < iterations_remainder ? 1 : 0);
        if (iterations_for_this_thread > 0) {
            threads.emplace_back(mccfr_task, iterations_for_this_thread, t);
        }
    }

    // Progress monitoring (simplified - no exploitability checks here)
    int report_interval_ms = 500;
    int iterations_reported_at_last_print = -1;
    bool first_report_triggered = false;

    if (num_iterations > 0) {
        std::cout << "MCCFR Iteration progress with exploitability tracking (target " << num_iterations << "):" << std::endl;
    }

    while(G_CFR_COMPLETED_ITERATIONS_COUNT.load(std::memory_order_acquire) < num_iterations) {
        std::this_thread::sleep_for(std::chrono::milliseconds(report_interval_ms));
        int current_completed = G_CFR_COMPLETED_ITERATIONS_COUNT.load(std::memory_order_acquire);
        
        if (!first_report_triggered || current_completed > iterations_reported_at_last_print) {
            std::cout << "  Completed " << current_completed << " / " << num_iterations << " iterations." << std::endl;
            iterations_reported_at_last_print = current_completed;
            first_report_triggered = true;
        }
        if (current_completed >= num_iterations) break; 
    }

    // Join threads
    for (auto& th : threads) {
        if (th.joinable()) {
            th.join();
        }
    }
    
    // Final completion message
    int final_completed_count = G_CFR_COMPLETED_ITERATIONS_COUNT.load(std::memory_order_acquire);
    if (num_iterations > 0) {
        if (final_completed_count > iterations_reported_at_last_print || !first_report_triggered) {
             std::cout << "  Completed " << final_completed_count << " / " << num_iterations << " iterations." << std::endl;
        }
        std::cout << "MCCFR Iteration " << final_completed_count << "/" << num_iterations << " (All threads complete)" << std::endl;
    }

    // Final exploitability computation
    std::cout << "\nComputing final exploitability..." << std::endl;
    try {
        double final_exploitability = compute_current_strategy_exploitability(start_state, exploitability_simulations * 2);
        double final_exploitability_percentage = (starting_pot > 0) ? (final_exploitability / starting_pot) * 100.0 : 0.0;
        {
            std::lock_guard<std::mutex> lock(exploitability_history_mutex);
            exploitability_history.push_back(std::make_tuple(final_completed_count, final_exploitability, final_exploitability_percentage));
        }
        std::cout << "Final exploitability: " << std::fixed << std::setprecision(6) << final_exploitability << " BB/hand" << std::endl;
        std::cout << "Final exploitability percentage: " << std::fixed << std::setprecision(2) << final_exploitability_percentage << "% of starting pot" << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Error computing final exploitability: " << e.what() << std::endl;
    }

    // Print exploitability history summary
    {
        std::lock_guard<std::mutex> lock(exploitability_history_mutex);
        if (!exploitability_history.empty()) {
            std::cout << "\nExploitability Training History:" << std::endl;
            std::cout << "Iteration\tExploitability (BB/hand)\tExploitability (%)" << std::endl;
            for (const auto& entry : exploitability_history) {
                std::cout << std::get<0>(entry) << "\t\t" << std::fixed << std::setprecision(6) << std::get<1>(entry) 
                          << "\t\t\t" << std::setprecision(2) << std::get<2>(entry) << "%" << std::endl;
            }
            
            if (exploitability_history.size() > 1) {
                double initial_exp = std::get<1>(exploitability_history[0]);
                double final_exp = std::get<1>(exploitability_history.back());
                double initial_exp_pct = std::get<2>(exploitability_history[0]);
                double final_exp_pct = std::get<2>(exploitability_history.back());
                double total_improvement = initial_exp - final_exp;
                double total_improvement_pct = initial_exp_pct - final_exp_pct;
                std::cout << "Total exploitability improvement: " << std::fixed << std::setprecision(6) 
                          << total_improvement << " BB/hand (" << std::fixed << std::setprecision(2) 
                          << (total_improvement / initial_exp * 100.0) << "% reduction)" << std::endl;
                std::cout << "Total exploitability percentage improvement: " << std::fixed << std::setprecision(2) 
                          << total_improvement_pct << " percentage points (" << std::fixed << std::setprecision(2) 
                          << (total_improvement_pct / initial_exp_pct * 100.0) << "% reduction)" << std::endl;
            }
        }
    }

    std::cout << "MCCFR iterations complete. Saving strategy..." << std::endl;
    std::ofstream outfile(output_csv_filename);
    outfile << "infoset_key,action_string,probability\n";

    for (const auto& pair : G_CFR_INFOSET_ACTIONS) {
        const std::string& infoset_key = pair.first;
        const auto& actions_for_infoset = pair.second;

        if (G_CFR_STRATEGY_SUM.count(infoset_key)) {
            const std::vector<double>& summed_strategy = G_CFR_STRATEGY_SUM.at(infoset_key);
            double total_summed_strategy = 0.0;
            for (double prob_sum : summed_strategy) {
                total_summed_strategy += prob_sum;
            }

            if (total_summed_strategy > 1e-9) {
                for (size_t a = 0; a < summed_strategy.size(); ++a) {
                    double avg_prob = summed_strategy[a] / total_summed_strategy;
                    outfile << "\"" << infoset_key << "\",\"" 
                            << action_to_string(actions_for_infoset[a].first, actions_for_infoset[a].second) 
                            << "\"," << std::fixed << std::setprecision(5) << avg_prob << "\n";
                }
            } else {
                if (!actions_for_infoset.empty()) {
                    double uniform_prob = 1.0 / static_cast<double>(actions_for_infoset.size());
                    for (size_t a = 0; a < actions_for_infoset.size(); ++a) {
                        outfile << "\"" << infoset_key << "\",\""
                                << action_to_string(actions_for_infoset[a].first, actions_for_infoset[a].second)
                                << "\"," << std::fixed << std::setprecision(5) << uniform_prob << "\n";
                    }
                }
            }
        } else {
            std::cerr << "Warning: Infoset " << infoset_key << " found in actions but not in strategy sum." << std::endl;
            if (!actions_for_infoset.empty()) {
                double uniform_prob = 1.0 / static_cast<double>(actions_for_infoset.size());
                 for (size_t a = 0; a < actions_for_infoset.size(); ++a) {
                    outfile << "\"" << infoset_key << "\",\""
                            << action_to_string(actions_for_infoset[a].first, actions_for_infoset[a].second)
                            << "\"," << std::fixed << std::setprecision(5) << uniform_prob << "\n";
                }
            }
        }
    }
    outfile.close();
    
    // Save exploitability history to file
    std::string exploitability_log_filename = output_csv_filename + "_exploitability_log.csv";
    std::ofstream exp_outfile(exploitability_log_filename);
    exp_outfile << "iteration,exploitability,exploitability_percentage\n";
    for (const auto& entry : exploitability_history) {
        exp_outfile << std::get<0>(entry) << "," << std::fixed << std::setprecision(6) << std::get<1>(entry) << "," << std::setprecision(2) << std::get<2>(entry) << "\n";
    }
    exp_outfile.close();
    
    std::cout << "Strategy saved to " << output_csv_filename << std::endl;
    std::cout << "Exploitability history saved to " << exploitability_log_filename << std::endl;
}

// Helper function to parse betting history string into actions
// Example: "c-r200-c" -> [CHECK, RAISE(200), CALL]
std::vector<std::pair<Action, double>> parse_betting_history(const std::string& betting_history) {
    std::vector<std::pair<Action, double>> actions;
    if (betting_history.empty()) return actions;
    
    std::istringstream iss(betting_history);
    std::string action_str;
    
    while (std::getline(iss, action_str, '-')) {
        if (action_str.empty()) continue;
        
        if (action_str == "c") {
            actions.push_back({Action::CALL, 0.0});
        } else if (action_str == "k" || action_str == "check") {
            actions.push_back({Action::CHECK, 0.0});
        } else if (action_str == "f" || action_str == "fold") {
            actions.push_back({Action::FOLD, 0.0});
        } else if (action_str[0] == 'r') {
            // Parse raise amount: "r200" -> RAISE with amount 200
            std::string amount_str = action_str.substr(1);
            double amount = std::stod(amount_str);
            actions.push_back({Action::BET_RAISE, amount});
        } else if (action_str[0] == 'b') {
            // Parse bet amount: "b100" -> BET with amount 100
            std::string amount_str = action_str.substr(1);
            double amount = std::stod(amount_str);
            actions.push_back({Action::BET_RAISE, amount});
        }
    }
    return actions;
}

// Helper function to parse community cards string
// Example: "As Kd Qh" -> [Card(A,s), Card(K,d), Card(Q,h)]
std::vector<Card> parse_community_cards(const std::string& cards_str) {
    std::vector<Card> cards;
    std::istringstream iss(cards_str);
    std::string card_str;
    
    while (iss >> card_str) {
        if (card_str.length() >= 2) {
            std::string rank = card_str.substr(0, card_str.length() - 1);
            char suit = card_str.back();
            cards.push_back(Card(rank, suit));
        }
    }
    return cards;
}

// Initialize a subgame state from betting history and community cards
HUState initialize_subgame_state(const HUGame& game_config, 
                                const std::string& betting_history,
                                const std::string& community_cards_str) {
    std::cout << "DEBUG: Initializing subgame state" << std::endl;
    std::cout << "DEBUG: Betting history: " << betting_history << std::endl;
    std::cout << "DEBUG: Community cards: " << community_cards_str << std::endl;

    HUState state(game_config);
    std::cout << "DEBUG: Created initial state" << std::endl;
    
    // Parse target community cards and remove them from deck before dealing
    std::vector<Card> target_community_cards = parse_community_cards(community_cards_str);
    std::cout << "DEBUG: Parsed " << target_community_cards.size() << " community cards" << std::endl;
    
    // Remove target community cards from deck before dealing hole cards
    if (!target_community_cards.empty()) {
        for (const auto& target_card : target_community_cards) {
            auto it = std::find_if(state.deck.begin(), state.deck.end(), 
                [&target_card](const Card& deck_card) {
                    return deck_card.rank == target_card.rank && deck_card.suit == target_card.suit;
                });
            if (it != state.deck.end()) {
                state.deck.erase(it);
                std::cout << "DEBUG: Removed " << target_card.toString() << " from deck" << std::endl;
            }
        }
        std::cout << "DEBUG: Deck size after removing community cards: " << state.deck.size() << std::endl;
    }
    
    state.apply_action(Action::DEAL); // Deal hole cards from modified deck
    std::cout << "DEBUG: Dealt hole cards" << std::endl;
    
    // Apply betting history step by step, letting the game naturally progress
    std::vector<std::pair<Action, double>> actions = parse_betting_history(betting_history);
    std::cout << "DEBUG: Parsed " << actions.size() << " betting actions" << std::endl;
    
    for (size_t i = 0; i < actions.size(); i++) {
        const auto& action_pair = actions[i];
        std::cout << "DEBUG: Processing action " << i + 1 << "/" << actions.size() 
                  << ", game_over=" << state.game_over << std::endl;
        
        if (state.game_over) {
            std::cout << "DEBUG: Game over, breaking action loop" << std::endl;
            break;
        }
        
        // Handle chance nodes (community card dealing) during natural game progression
        while (state.is_chance_node() && !state.game_over) {
            std::cout << "DEBUG: Handling chance node, round=" << state.round << std::endl;
            state.apply_action(Action::DEAL);
        }
        
        // Apply the betting action if game is not over and not a chance node
        if (!state.game_over && !state.is_chance_node()) {
            std::cout << "DEBUG: Applying betting action, round=" << state.round << std::endl;
            state.apply_action(action_pair.first, action_pair.second);
        }
    }
    
    // Now set the community cards to match the target
    if (!target_community_cards.empty()) {
        std::cout << "DEBUG: Setting target community cards" << std::endl;
        state.community_cards = target_community_cards;
        
        // Update the round based on the number of target community cards
        if (target_community_cards.size() == 0) {
            state.round = "preflop";
        } else if (target_community_cards.size() == 3) {
            state.round = "flop";
        } else if (target_community_cards.size() == 4) {
            state.round = "turn";
        } else if (target_community_cards.size() == 5) {
            state.round = "river";
        }
        std::cout << "DEBUG: Updated round to " << state.round << std::endl;
    }
    
    std::cout << "DEBUG: Final state - round=" << state.round 
              << ", pot=" << state.pot 
              << ", community_cards=" << state.community_cards.size() << std::endl;
    
    return state;
}

// Get all possible hole card combinations for a player, excluding known cards
std::vector<std::pair<Card, Card>> get_possible_hole_cards(const std::vector<Card>& known_cards) {
    std::vector<std::pair<Card, Card>> possible_hands;
    std::vector<std::string> ranks = {"2", "3", "4", "5", "6", "7", "8", "9", "10", "J", "Q", "K", "A"};
    std::vector<std::string> suits = {"c", "d", "h", "s"};
    
    // Create deck excluding known cards
    std::vector<Card> available_cards;
    for (const auto& rank : ranks) {
        for (const auto& suit : suits) {
            Card card(rank, suit[0]);
            bool is_known = false;
            // Check if card is in known_cards
            for (const auto& known_card : known_cards) {
                if (card.rank == known_card.rank && card.suit == known_card.suit) {
                    is_known = true;
                    break;
                }
            }
            if (!is_known) {
                available_cards.push_back(card);
            }
        }
    }
    
    // Generate all possible 2-card combinations
    for (size_t i = 0; i < available_cards.size(); ++i) {
        for (size_t j = i + 1; j < available_cards.size(); ++j) {
            possible_hands.push_back({available_cards[i], available_cards[j]});
        }
    }
    
    return possible_hands;
}

// Sample remaining community cards for chance nodes
std::vector<Card> sample_remaining_community_cards(const std::vector<Card>& current_community,
                                                  const std::vector<Card>& known_cards,
                                                  int target_community_size,
                                                  std::mt19937& rng) {
    std::vector<Card> result = current_community;
    if (result.size() >= static_cast<size_t>(target_community_size)) {
        return result;
    }
    
    // Create deck excluding known cards
    std::vector<std::string> ranks = {"2", "3", "4", "5", "6", "7", "8", "9", "10", "J", "Q", "K", "A"};
    std::vector<std::string> suits = {"c", "d", "h", "s"};
    std::vector<Card> available_cards;
    
    for (const auto& rank : ranks) {
        for (const auto& suit : suits) {
            Card card(rank, suit[0]);
            bool is_known = false;
            for (const auto& known_card : known_cards) {
                if (card.rank == known_card.rank && card.suit == known_card.suit) {
                    is_known = true;
                    break;
                }
            }
            if (!is_known) {
                available_cards.push_back(card);
            }
        }
    }
    
    // Sample needed cards
    std::shuffle(available_cards.begin(), available_cards.end(), rng);
    int cards_needed = target_community_size - result.size();
    for (int i = 0; i < cards_needed && i < static_cast<int>(available_cards.size()); ++i) {
        result.push_back(available_cards[i]);
    }
    
    return result;
}

// Modified CFR Recursive Function for Subgame Solving
std::vector<double> cfr_subgame_recursive(HUState current_state, 
                                         std::vector<double> reach_probs,
                                         const std::vector<std::pair<Card, Card>>& possible_hands_p0,
                                         const std::vector<std::pair<Card, Card>>& possible_hands_p1,
                                         int sampled_hand_idx_p0,
                                         int sampled_hand_idx_p1) {
    
    // DEBUG: Add entry debug
    static std::atomic<int> call_count(0);
    int this_call = call_count.fetch_add(1);
    if (this_call < 5) { // Only log first few calls to avoid spam
        std::cout << "DEBUG: cfr_subgame_recursive call #" << this_call 
                  << ", game_over=" << current_state.game_over 
                  << ", is_chance=" << current_state.is_chance_node()
                  << ", round=" << current_state.round << std::endl;
    }
    
    if (current_state.game_over) {
        if (this_call < 5) {
            std::cout << "DEBUG: Game over, returning: " << current_state.returns()[0] << ", " << current_state.returns()[1] << std::endl;
        }
        return current_state.returns();
    }

    if (current_state.is_chance_node()) {
        if (this_call < 5) {
            std::cout << "DEBUG: Handling chance node" << std::endl;
        }
        // Handle chance nodes by sampling remaining community cards if needed
        HUState next_state = current_state;
        
        // Check if we need more community cards
        int target_community_size = 0;
        if (current_state.round == "flop" && current_state.community_cards.size() < 3) {
            target_community_size = 3;
        } else if (current_state.round == "turn" && current_state.community_cards.size() < 4) {
            target_community_size = 4;
        } else if (current_state.round == "river" && current_state.community_cards.size() < 5) {
            target_community_size = 5;
        }
        
        if (target_community_size > 0) {
            // Exclude both the original community cards AND the sampled hole cards for both players
            std::vector<Card> cards_to_exclude = G_ORIGINAL_COMMUNITY_CARDS;
            
            // Add sampled hole cards to exclusion list - check bounds first
            if (sampled_hand_idx_p0 >= 0 && sampled_hand_idx_p0 < static_cast<int>(possible_hands_p0.size())) {
                cards_to_exclude.push_back(possible_hands_p0[sampled_hand_idx_p0].first);
                cards_to_exclude.push_back(possible_hands_p0[sampled_hand_idx_p0].second);
            }
            if (sampled_hand_idx_p1 >= 0 && sampled_hand_idx_p1 < static_cast<int>(possible_hands_p1.size())) {
                cards_to_exclude.push_back(possible_hands_p1[sampled_hand_idx_p1].first);
                cards_to_exclude.push_back(possible_hands_p1[sampled_hand_idx_p1].second);
            }
            
            std::lock_guard<std::mutex> rng_lock(G_SUBGAME_RNG_MUTEX);
            next_state.community_cards = sample_remaining_community_cards(
                current_state.community_cards, cards_to_exclude, target_community_size, G_SUBGAME_RNG);
        }
        
        next_state.apply_action(Action::DEAL);
        return cfr_subgame_recursive(next_state, reach_probs, possible_hands_p0, possible_hands_p1, 
                                   sampled_hand_idx_p0, sampled_hand_idx_p1);
    }

    // DEBUG: Player decision node
    if (this_call < 5) {
        std::cout << "DEBUG: Player decision node, player=" << current_state.current_player() << std::endl;
    }

    // Set the sampled hole cards in the state for infoset calculation - check bounds first
    HUState state_with_cards = current_state;
    if (sampled_hand_idx_p0 >= 0 && sampled_hand_idx_p0 < static_cast<int>(possible_hands_p0.size())) {
        if (state_with_cards.cards.size() >= 2) {
            state_with_cards.cards[0] = possible_hands_p0[sampled_hand_idx_p0].first;
            state_with_cards.cards[1] = possible_hands_p0[sampled_hand_idx_p0].second;
        }
    }
    if (sampled_hand_idx_p1 >= 0 && sampled_hand_idx_p1 < static_cast<int>(possible_hands_p1.size())) {
        if (state_with_cards.cards.size() >= 4) {
            state_with_cards.cards[2] = possible_hands_p1[sampled_hand_idx_p1].first;
            state_with_cards.cards[3] = possible_hands_p1[sampled_hand_idx_p1].second;
        }
    }

    int player_idx = current_state.current_player();
    std::string infoset_key = get_infoset_key_for_state(state_with_cards);

    // DEBUG: Infoset key
    if (this_call < 5) {
        std::cout << "DEBUG: Generated infoset key: '" << infoset_key << "'" << std::endl;
    }

    std::vector<std::pair<Action, double>> legal_actions;
    int num_actions_for_infoset;

    { // Scope for G_CFR_INFOSET_ACTIONS_MUTEX
        std::lock_guard<std::mutex> lock(G_CFR_INFOSET_ACTIONS_MUTEX);
        auto it_actions = G_CFR_INFOSET_ACTIONS.find(infoset_key);
        if (it_actions == G_CFR_INFOSET_ACTIONS.end()) {
            legal_actions = current_state.legal_actions(
                current_state.game.pot_fraction_bet_sizes,
                current_state.game.fixed_bet_sizes_bb
            );
            
            // DEBUG: Legal actions
            if (this_call < 5) {
                std::cout << "DEBUG: Found " << legal_actions.size() << " legal actions for new infoset" << std::endl;
            }
            
            if (legal_actions.empty()) {
                std::cerr << "CFR Warning: Player " << player_idx << " has no legal actions at infoset "
                          << infoset_key << " but game not over and not chance node. State:\n"
                          << current_state.to_string(false) << std::endl;
                return current_state.returns();
            }
            G_CFR_INFOSET_ACTIONS[infoset_key] = legal_actions;
            G_CFR_REGRET_SUM.try_emplace(infoset_key, std::vector<double>(legal_actions.size(), 0.0));
            G_CFR_STRATEGY_SUM.try_emplace(infoset_key, std::vector<double>(legal_actions.size(), 0.0));
            num_actions_for_infoset = legal_actions.size();
            
            // DEBUG: New infoset registered
            if (this_call < 5) {
                std::cout << "DEBUG: Registered new infoset with " << num_actions_for_infoset << " actions" << std::endl;
            }
        } else {
            legal_actions = it_actions->second;
            num_actions_for_infoset = legal_actions.size();
            
            // DEBUG: Existing infoset
            if (this_call < 5) {
                std::cout << "DEBUG: Found existing infoset with " << num_actions_for_infoset << " actions" << std::endl;
            }
        }
    } // G_CFR_INFOSET_ACTIONS_MUTEX is released

    if (num_actions_for_infoset == 0) { 
         std::cerr << "CFR Error: num_actions is 0 for infoset " << infoset_key << std::endl;
         return current_state.returns();
    }

    std::vector<double> strategy = get_mccfr_strategy(infoset_key, num_actions_for_infoset);
    std::vector<double> action_utilities(num_actions_for_infoset);
    std::vector<double> node_expected_value(2, 0.0); // EV for P0, P1 from this node

    for (int a = 0; a < num_actions_for_infoset; ++a) {
        HUState next_state = current_state;
        next_state.apply_action(legal_actions[a].first, legal_actions[a].second);
        
        std::vector<double> child_reach_probs = reach_probs;
        child_reach_probs[player_idx] *= strategy[a];
        
        std::vector<double> child_node_utilities = cfr_subgame_recursive(next_state, child_reach_probs, possible_hands_p0, possible_hands_p1, 
                                   sampled_hand_idx_p0, sampled_hand_idx_p1);
        action_utilities[a] = child_node_utilities[player_idx];
        
        for(int p=0; p<2; ++p) {
            node_expected_value[p] += strategy[a] * child_node_utilities[p];
        }
    }

    // Update regrets for the current player (player_idx)
    double opponent_reach_prob = reach_probs[1 - player_idx];
    { // Scope for G_CFR_REGRET_SUM_MUTEX
        std::lock_guard<std::mutex> lock(G_CFR_REGRET_SUM_MUTEX);
        
        // Ensure infoset exists in regret sum
        if (G_CFR_REGRET_SUM.find(infoset_key) == G_CFR_REGRET_SUM.end()) {
            G_CFR_REGRET_SUM[infoset_key] = std::vector<double>(num_actions_for_infoset, 0.0);
        }
        
        // Ensure correct size
        if (static_cast<int>(G_CFR_REGRET_SUM[infoset_key].size()) != num_actions_for_infoset) {
            G_CFR_REGRET_SUM[infoset_key].resize(num_actions_for_infoset, 0.0);
        }
        
        for (int a = 0; a < num_actions_for_infoset; ++a) {
            double regret = action_utilities[a] - node_expected_value[player_idx];
            G_CFR_REGRET_SUM[infoset_key][a] += opponent_reach_prob * regret;
        }
    } // G_CFR_REGRET_SUM_MUTEX is released
    
    // Update average strategy sum for the current player (player_idx)
    double self_reach_prob = reach_probs[player_idx];
    { // Scope for G_CFR_STRATEGY_SUM_MUTEX
        std::lock_guard<std::mutex> lock(G_CFR_STRATEGY_SUM_MUTEX);
        
        // Ensure infoset exists in strategy sum
        if (G_CFR_STRATEGY_SUM.find(infoset_key) == G_CFR_STRATEGY_SUM.end()) {
            G_CFR_STRATEGY_SUM[infoset_key] = std::vector<double>(num_actions_for_infoset, 0.0);
        }
        
        // Ensure correct size
        if (static_cast<int>(G_CFR_STRATEGY_SUM[infoset_key].size()) != num_actions_for_infoset) {
            G_CFR_STRATEGY_SUM[infoset_key].resize(num_actions_for_infoset, 0.0);
        }
        
        // Linear CFR: Weight strategy updates by iteration number for faster convergence
        double weight = G_USE_LINEAR_CFR ? std::max(1, G_CFR_CURRENT_ITERATION.load()) : 1.0;
        
        for (int a = 0; a < num_actions_for_infoset; ++a) {
            G_CFR_STRATEGY_SUM[infoset_key][a] += weight * strategy[a];
        }
    } // G_CFR_STRATEGY_SUM_MUTEX is released
    
    return node_expected_value;
}

// Helper function to get possible hole cards for a player, excluding specific cards
std::vector<std::pair<Card, Card>> get_possible_hole_cards_excluding(const std::vector<Card>& cards_to_exclude) {
    std::vector<std::pair<Card, Card>> possible_hands;
    std::vector<std::string> ranks = {"2", "3", "4", "5", "6", "7", "8", "9", "10", "J", "Q", "K", "A"};
    std::vector<std::string> suits = {"c", "d", "h", "s"};
    
    // Create deck excluding specified cards
    std::vector<Card> available_cards;
    for (const auto& rank : ranks) {
        for (const auto& suit : suits) {
            Card card(rank, suit[0]);
            bool is_excluded = false;
            // Check if card is in cards_to_exclude
            for (const auto& excluded_card : cards_to_exclude) {
                if (card.rank == excluded_card.rank && card.suit == excluded_card.suit) {
                    is_excluded = true;
                    break;
                }
            }
            if (!is_excluded) {
                available_cards.push_back(card);
            }
        }
    }
    
    // Generate all possible 2-card combinations
    for (size_t i = 0; i < available_cards.size(); ++i) {
        for (size_t j = i + 1; j < available_cards.size(); ++j) {
            possible_hands.push_back({available_cards[i], available_cards[j]});
        }
    }
    
    return possible_hands;
}

// Helper function to check if two hands overlap (share any cards)
bool hands_overlap(const std::pair<Card, Card>& hand1, const std::pair<Card, Card>& hand2) {
    return (hand1.first.rank == hand2.first.rank && hand1.first.suit == hand2.first.suit) ||
           (hand1.first.rank == hand2.second.rank && hand1.first.suit == hand2.second.suit) ||
           (hand1.second.rank == hand2.first.rank && hand1.second.suit == hand2.first.suit) ||
           (hand1.second.rank == hand2.second.rank && hand1.second.suit == hand2.second.suit);
}

// Main function to solve subgame using CFR
void solve_subgame_cfr(const HUGame& game_config,
                      const std::string& betting_history,
                      const std::string& community_cards_str,
                      int num_iterations,
                      const std::string& output_csv_filename,
                      int num_threads_to_use = 1) {
    
    // Clear global CFR data structures
    G_CFR_REGRET_SUM.clear();
    G_CFR_STRATEGY_SUM.clear();
    G_CFR_INFOSET_ACTIONS.clear();
    G_CFR_COMPLETED_ITERATIONS_COUNT = 0;

    std::cout << "Initializing subgame state..." << std::endl;
    std::cout << "Betting history: " << betting_history << std::endl;
    std::cout << "Community cards: " << community_cards_str << std::endl;

    HUState subgame_root = initialize_subgame_state(game_config, betting_history, community_cards_str);
    std::cout << "Subgame root state:\n" << subgame_root.to_string(false) << std::endl;
    
    // DEBUG: Check subgame root state details
    std::cout << "DEBUG: Subgame root - game_over=" << subgame_root.game_over 
              << ", is_chance=" << subgame_root.is_chance_node()
              << ", current_player=" << subgame_root.current_player()
              << ", round=" << subgame_root.round << std::endl;

    // Store the original community cards that were explicitly input by the user
    G_ORIGINAL_COMMUNITY_CARDS = subgame_root.community_cards;

    // Get possible hole card combinations for both players (excluding community cards)
    std::vector<Card> known_cards = subgame_root.community_cards;
    std::vector<std::pair<Card, Card>> possible_hands_p0 = get_possible_hole_cards_excluding(known_cards);
    std::vector<std::pair<Card, Card>> possible_hands_p1 = get_possible_hole_cards_excluding(known_cards);
    
    std::cout << "Found " << possible_hands_p0.size() << " possible hands for P0" << std::endl;
    std::cout << "Found " << possible_hands_p1.size() << " possible hands for P1" << std::endl;
    std::cout << "Starting subgame CFR for " << num_iterations << " iterations..." << std::endl;

    std::vector<std::thread> threads;
    int iterations_per_thread_base = num_iterations / num_threads_to_use;
    int iterations_remainder = num_iterations % num_threads_to_use;

    auto cfr_subgame_task = [&](int iterations_for_this_thread, int thread_id) {
        std::mt19937 thread_rng(std::random_device{}() + thread_id);
        
        // DEBUG: Task start
        std::cout << "DEBUG: Thread " << thread_id << " starting with " << iterations_for_this_thread << " iterations" << std::endl;
        
        for (int i = 0; i < iterations_for_this_thread; ++i) {
            // DEBUG: Iteration start
            if (i < 3) { // Only log first few iterations per thread
                std::cout << "DEBUG: Thread " << thread_id << " iteration " << i << std::endl;
            }
            
            // Sample hole cards for P0
            std::uniform_int_distribution<int> hand_dist_p0(0, possible_hands_p0.size() - 1);
            int sampled_hand_idx_p0 = hand_dist_p0(thread_rng);
            
            // DEBUG: Sampled P0 hand
            if (i < 3) {
                std::cout << "DEBUG: Thread " << thread_id << " sampled P0 hand idx=" << sampled_hand_idx_p0 
                          << " (" << possible_hands_p0[sampled_hand_idx_p0].first.toString() 
                          << " " << possible_hands_p0[sampled_hand_idx_p0].second.toString() << ")" << std::endl;
            }
            
            // Sample hole cards for P1, ensuring no overlap with P0's cards
            std::uniform_int_distribution<int> hand_dist_p1(0, possible_hands_p1.size() - 1);
            int sampled_hand_idx_p1;
            int max_attempts = 100; // Prevent infinite loop
            int attempts = 0;
            
            do {
                sampled_hand_idx_p1 = hand_dist_p1(thread_rng);
                attempts++;
            } while (attempts < max_attempts && 
                     hands_overlap(possible_hands_p0[sampled_hand_idx_p0], possible_hands_p1[sampled_hand_idx_p1]));
            
            if (attempts >= max_attempts) {
                std::cerr << "Warning: Could not find non-overlapping hands after " << max_attempts << " attempts" << std::endl;
                continue; // Skip this iteration
            }
            
            // DEBUG: Sampled P1 hand
            if (i < 3) {
                std::cout << "DEBUG: Thread " << thread_id << " sampled P1 hand idx=" << sampled_hand_idx_p1 
                          << " (" << possible_hands_p1[sampled_hand_idx_p1].first.toString() 
                          << " " << possible_hands_p1[sampled_hand_idx_p1].second.toString() << ")" << std::endl;
            }
            
            // DEBUG: About to call CFR recursive
            if (i < 3) {
                std::cout << "DEBUG: Thread " << thread_id << " calling cfr_subgame_recursive for iteration " << i << std::endl;
            }
            
            std::vector<double> initial_reach_probs(2, 1.0);
            cfr_subgame_recursive(subgame_root, initial_reach_probs, 
                                possible_hands_p0, possible_hands_p1,
                                sampled_hand_idx_p0, sampled_hand_idx_p1);
            
            // DEBUG: CFR recursive completed
            if (i < 3) {
                std::cout << "DEBUG: Thread " << thread_id << " completed cfr_subgame_recursive for iteration " << i << std::endl;
            }
            
            G_CFR_COMPLETED_ITERATIONS_COUNT.fetch_add(1, std::memory_order_relaxed);
        }
        
        // DEBUG: Task end
        std::cout << "DEBUG: Thread " << thread_id << " completed all iterations" << std::endl;
    };

    for (int t = 0; t < num_threads_to_use; ++t) {
        int iterations_for_this_thread = iterations_per_thread_base + (t < iterations_remainder ? 1 : 0);
        if (iterations_for_this_thread > 0) {
            threads.emplace_back(cfr_subgame_task, iterations_for_this_thread, t);
        }
    }

    // Progress monitoring
    int report_interval_ms = 500;
    int iterations_reported_at_last_print = -1;
    bool first_report_triggered = false;

    if (num_iterations > 0) {
        std::cout << "Subgame CFR iteration progress (target " << num_iterations << "):" << std::endl;
    }

    while(G_CFR_COMPLETED_ITERATIONS_COUNT.load(std::memory_order_acquire) < num_iterations) {
        std::this_thread::sleep_for(std::chrono::milliseconds(report_interval_ms));
        int current_completed = G_CFR_COMPLETED_ITERATIONS_COUNT.load(std::memory_order_acquire);
        
        if (!first_report_triggered || current_completed > iterations_reported_at_last_print) {
            std::cout << "  Completed " << current_completed << " / " << num_iterations << " iterations." << std::endl;
            iterations_reported_at_last_print = current_completed;
            first_report_triggered = true;
        }
        if (current_completed >= num_iterations) break; 
    }

    // Join threads
    for (auto& th : threads) {
        if (th.joinable()) {
            th.join();
        }
    }
    
    std::cout << "Subgame CFR iterations complete. Saving strategy..." << std::endl;
    
    // DEBUG: Check how many infosets were found
    std::cout << "DEBUG: Found " << G_CFR_INFOSET_ACTIONS.size() << " infosets total" << std::endl;
    std::cout << "DEBUG: G_CFR_STRATEGY_SUM size: " << G_CFR_STRATEGY_SUM.size() << std::endl;
    std::cout << "DEBUG: G_CFR_REGRET_SUM size: " << G_CFR_REGRET_SUM.size() << std::endl;
    
    // Save strategy to CSV
    std::ofstream outfile(output_csv_filename);
    outfile << "infoset_key,action_string,probability\n";

    int lines_written = 0;
    for (const auto& pair : G_CFR_INFOSET_ACTIONS) {
        const std::string& infoset_key = pair.first;
        const auto& actions_for_infoset = pair.second;

        // DEBUG: Processing infoset
        if (lines_written < 5) {
            std::cout << "DEBUG: Processing infoset '" << infoset_key << "' with " << actions_for_infoset.size() << " actions" << std::endl;
        }

        if (G_CFR_STRATEGY_SUM.count(infoset_key)) {
            const std::vector<double>& summed_strategy = G_CFR_STRATEGY_SUM.at(infoset_key);
            double total_summed_strategy = 0.0;
            for (double prob_sum : summed_strategy) {
                total_summed_strategy += prob_sum;
            }
            
            // DEBUG: Strategy sum details
            if (lines_written < 5) {
                std::cout << "DEBUG: Infoset '" << infoset_key << "' total_summed_strategy=" << total_summed_strategy << std::endl;
            }

            if (total_summed_strategy > 1e-9) {
                for (size_t a = 0; a < summed_strategy.size(); ++a) {
                    double avg_prob = summed_strategy[a] / total_summed_strategy;
                    outfile << "\"" << infoset_key << "\",\"" 
                            << action_to_string(actions_for_infoset[a].first, actions_for_infoset[a].second) 
                            << "\"," << std::fixed << std::setprecision(5) << avg_prob << "\n";
                    lines_written++;
                }
            } else {
                if (!actions_for_infoset.empty()) {
                    double uniform_prob = 1.0 / static_cast<double>(actions_for_infoset.size());
                    for (size_t a = 0; a < actions_for_infoset.size(); ++a) {
                        outfile << "\"" << infoset_key << "\",\""
                                << action_to_string(actions_for_infoset[a].first, actions_for_infoset[a].second)
                                << "\"," << std::fixed << std::setprecision(5) << uniform_prob << "\n";
                        lines_written++;
                    }
                }
            }
        }
    }
    outfile.close();
    
    // DEBUG: Final summary
    std::cout << "DEBUG: Wrote " << lines_written << " lines to CSV file" << std::endl;
    std::cout << "Subgame strategy saved to " << output_csv_filename << std::endl;
}

// Helper function to convert action history from HUState to betting history string
std::string extract_betting_history_from_state(const HUState& state) {
    std::vector<std::string> action_strings;
    
    // Iterate through all rounds in the action history
    std::vector<std::string> round_order = {"preflop", "flop", "turn", "river"};
    
    for (const std::string& round_name : round_order) {
        auto it = state.round_action_history.find(round_name);
        if (it != state.round_action_history.end()) {
            const auto& round_actions = it->second;
            for (const auto& action_tuple : round_actions) {
                Action action = std::get<1>(action_tuple);
                double amount = std::get<2>(action_tuple);
                
                // Convert action to betting history format
                std::string action_str;
                switch (action) {
                    case Action::CHECK:
                        action_str = "k";
                        break;
                    case Action::CALL:
                        action_str = "c";
                        break;
                    case Action::FOLD:
                        action_str = "f";
                        break;
                    case Action::BET_RAISE:
                    case Action::ALL_IN:
                        action_str = "r" + std::to_string(static_cast<int>(std::round(amount)));
                        break;
                    case Action::POST_SB:
                    case Action::POST_BB:
                        // Skip blind posts as they're not part of betting history
                        continue;
                    default:
                        continue;
                }
                action_strings.push_back(action_str);
            }
        }
        
        // Only include actions up to current round
        if (round_name == state.round) break;
    }
    
    // Join actions with "-" separator
    std::string result;
    for (size_t i = 0; i < action_strings.size(); ++i) {
        if (i > 0) result += "-";
        result += action_strings[i];
    }
    return result;
}

// Modified CFR Recursive Function for Subgame Solving using External Sampling
double mccfr_subgame_recursive(HUState current_state, 
                              int updating_player,
                              const std::vector<std::pair<Card, Card>>& possible_hands_p0,
                              const std::vector<std::pair<Card, Card>>& possible_hands_p1,
                              int sampled_hand_idx_p0,
                              int sampled_hand_idx_p1,
                              std::mt19937& rng) {
    
    static std::atomic<int> call_count(0);
    int this_call = call_count.fetch_add(1);
    if (this_call < 5) {
        std::cout << "DEBUG: mccfr_subgame_recursive call #" << this_call 
                  << ", game_over=" << current_state.game_over 
                  << ", is_chance=" << current_state.is_chance_node()
                  << ", round=" << current_state.round << std::endl;
    }
    
    if (current_state.game_over) {
        if (this_call < 5) {
            std::cout << "DEBUG: Game over, returning utility for player " << updating_player << std::endl;
        }
        std::vector<double> returns = current_state.returns();
        return returns[updating_player];
    }

    if (current_state.is_chance_node()) {
        if (this_call < 5) {
            std::cout << "DEBUG: Handling chance node" << std::endl;
        }
        HUState next_state = current_state;
        
        // Check if we need more community cards
        int target_community_size = 0;
        if (current_state.round == "flop" && current_state.community_cards.size() < 3) {
            target_community_size = 3;
        } else if (current_state.round == "turn" && current_state.community_cards.size() < 4) {
            target_community_size = 4;
        } else if (current_state.round == "river" && current_state.community_cards.size() < 5) {
            target_community_size = 5;
        }
        
        if (target_community_size > 0) {
            std::vector<Card> cards_to_exclude = G_ORIGINAL_COMMUNITY_CARDS;
            
            if (sampled_hand_idx_p0 >= 0 && sampled_hand_idx_p0 < static_cast<int>(possible_hands_p0.size())) {
                cards_to_exclude.push_back(possible_hands_p0[sampled_hand_idx_p0].first);
                cards_to_exclude.push_back(possible_hands_p0[sampled_hand_idx_p0].second);
            }
            if (sampled_hand_idx_p1 >= 0 && sampled_hand_idx_p1 < static_cast<int>(possible_hands_p1.size())) {
                cards_to_exclude.push_back(possible_hands_p1[sampled_hand_idx_p1].first);
                cards_to_exclude.push_back(possible_hands_p1[sampled_hand_idx_p1].second);
            }
            
            std::lock_guard<std::mutex> rng_lock(G_SUBGAME_RNG_MUTEX);
            next_state.community_cards = sample_remaining_community_cards(
                current_state.community_cards, cards_to_exclude, target_community_size, G_SUBGAME_RNG);
        }
        
        next_state.apply_action(Action::DEAL);
        return mccfr_subgame_recursive(next_state, updating_player, possible_hands_p0, possible_hands_p1, 
                                      sampled_hand_idx_p0, sampled_hand_idx_p1, rng);
    }

    if (this_call < 5) {
        std::cout << "DEBUG: Player decision node, player=" << current_state.current_player() << std::endl;
    }

    // Set the sampled hole cards in the state for infoset calculation
    HUState state_with_cards = current_state;
    if (sampled_hand_idx_p0 >= 0 && sampled_hand_idx_p0 < static_cast<int>(possible_hands_p0.size())) {
        if (state_with_cards.cards.size() >= 2) {
            state_with_cards.cards[0] = possible_hands_p0[sampled_hand_idx_p0].first;
            state_with_cards.cards[1] = possible_hands_p0[sampled_hand_idx_p0].second;
        }
    }
    if (sampled_hand_idx_p1 >= 0 && sampled_hand_idx_p1 < static_cast<int>(possible_hands_p1.size())) {
        if (state_with_cards.cards.size() >= 4) {
            state_with_cards.cards[2] = possible_hands_p1[sampled_hand_idx_p1].first;
            state_with_cards.cards[3] = possible_hands_p1[sampled_hand_idx_p1].second;
        }
    }

    int current_player = current_state.current_player();
    std::string infoset_key = get_infoset_key_for_state(state_with_cards);

    if (this_call < 5) {
        std::cout << "DEBUG: Generated infoset key: '" << infoset_key << "'" << std::endl;
    }

    std::vector<std::pair<Action, double>> legal_actions;
    int num_actions_for_infoset;

    { // Scope for G_CFR_INFOSET_ACTIONS_MUTEX
        std::lock_guard<std::mutex> lock(G_CFR_INFOSET_ACTIONS_MUTEX);
        auto it_actions = G_CFR_INFOSET_ACTIONS.find(infoset_key);
        if (it_actions == G_CFR_INFOSET_ACTIONS.end()) {
            legal_actions = current_state.legal_actions(
                current_state.game.pot_fraction_bet_sizes,
                current_state.game.fixed_bet_sizes_bb
            );
            
            if (this_call < 5) {
                std::cout << "DEBUG: Found " << legal_actions.size() << " legal actions for new infoset" << std::endl;
            }
            
            if (legal_actions.empty()) {
                std::cerr << "MCCFR Warning: Player " << current_player << " has no legal actions at infoset "
                          << infoset_key << " but game not over and not chance node." << std::endl;
                std::vector<double> returns = current_state.returns();
                return returns[updating_player];
            }
            G_CFR_INFOSET_ACTIONS[infoset_key] = legal_actions;
            G_CFR_REGRET_SUM.try_emplace(infoset_key, std::vector<double>(legal_actions.size(), 0.0));
            G_CFR_STRATEGY_SUM.try_emplace(infoset_key, std::vector<double>(legal_actions.size(), 0.0));
            num_actions_for_infoset = legal_actions.size();
            
            if (this_call < 5) {
                std::cout << "DEBUG: Registered new infoset with " << num_actions_for_infoset << " actions" << std::endl;
            }
        } else {
            legal_actions = it_actions->second;
            num_actions_for_infoset = legal_actions.size();
            
            if (this_call < 5) {
                std::cout << "DEBUG: Found existing infoset with " << num_actions_for_infoset << " actions" << std::endl;
            }
        }
    }

    if (num_actions_for_infoset == 0) { 
         std::cerr << "MCCFR Error: num_actions is 0 for infoset " << infoset_key << std::endl;
         std::vector<double> returns = current_state.returns();
         return returns[updating_player];
    }

    std::vector<double> strategy = get_mccfr_strategy(infoset_key, num_actions_for_infoset);

    if (current_player == updating_player) {
        // Player being updated: traverse all actions
        std::vector<double> action_utilities(num_actions_for_infoset);
        double node_utility = 0.0;

        for (int a = 0; a < num_actions_for_infoset; ++a) {
            HUState next_state = current_state;
            next_state.apply_action(legal_actions[a].first, legal_actions[a].second);
            
            action_utilities[a] = mccfr_subgame_recursive(
                next_state, updating_player, possible_hands_p0, possible_hands_p1, 
                sampled_hand_idx_p0, sampled_hand_idx_p1, rng);
            node_utility += strategy[a] * action_utilities[a];
        }

        // Update regrets
        { 
            std::lock_guard<std::mutex> lock(G_CFR_REGRET_SUM_MUTEX);
            
            // Ensure infoset exists in regret sum
            if (G_CFR_REGRET_SUM.find(infoset_key) == G_CFR_REGRET_SUM.end()) {
                G_CFR_REGRET_SUM[infoset_key] = std::vector<double>(num_actions_for_infoset, 0.0);
            }
            
            // Ensure correct size
            if (static_cast<int>(G_CFR_REGRET_SUM[infoset_key].size()) != num_actions_for_infoset) {
                G_CFR_REGRET_SUM[infoset_key].resize(num_actions_for_infoset, 0.0);
            }
            
            for (int a = 0; a < num_actions_for_infoset; ++a) {
                double regret = action_utilities[a] - node_utility;
                G_CFR_REGRET_SUM[infoset_key][a] += regret;
            }
        }
        
        // Update average strategy
        { 
            std::lock_guard<std::mutex> lock(G_CFR_STRATEGY_SUM_MUTEX);
            
            // Ensure infoset exists in strategy sum
            if (G_CFR_STRATEGY_SUM.find(infoset_key) == G_CFR_STRATEGY_SUM.end()) {
                G_CFR_STRATEGY_SUM[infoset_key] = std::vector<double>(num_actions_for_infoset, 0.0);
            }
            
            // Ensure correct size
            if (static_cast<int>(G_CFR_STRATEGY_SUM[infoset_key].size()) != num_actions_for_infoset) {
                G_CFR_STRATEGY_SUM[infoset_key].resize(num_actions_for_infoset, 0.0);
            }
            
            // Linear CFR: Weight strategy updates by iteration number for faster convergence
            double weight = G_USE_LINEAR_CFR ? std::max(1, G_CFR_CURRENT_ITERATION.load()) : 1.0;
            
            for (int a = 0; a < num_actions_for_infoset; ++a) {
                G_CFR_STRATEGY_SUM[infoset_key][a] += weight * strategy[a];
            }
        }
        
        return node_utility;
    } else {
        // Opponent: sample action according to strategy
        std::discrete_distribution<int> action_dist(strategy.begin(), strategy.end());
        int sampled_action = action_dist(rng);
        
        HUState next_state = current_state;
        next_state.apply_action(legal_actions[sampled_action].first, legal_actions[sampled_action].second);
        
        return mccfr_subgame_recursive(
            next_state, updating_player, possible_hands_p0, possible_hands_p1, 
            sampled_hand_idx_p0, sampled_hand_idx_p1, rng);
    }
}

// Main function to solve subgame using MCCFR
void solve_subgame_mccfr(const HUGame& game_config,
                        const std::string& betting_history,
                        const std::string& community_cards_str,
                        int num_iterations,
                        const std::string& output_csv_filename,
                        int num_threads_to_use = 1) {
    
    // Clear global CFR data structures
    G_CFR_REGRET_SUM.clear();
    G_CFR_STRATEGY_SUM.clear();
    G_CFR_INFOSET_ACTIONS.clear();
    G_CFR_COMPLETED_ITERATIONS_COUNT = 0;

    std::cout << "Initializing subgame state..." << std::endl;
    std::cout << "Betting history: " << betting_history << std::endl;
    std::cout << "Community cards: " << community_cards_str << std::endl;

    HUState subgame_root = initialize_subgame_state(game_config, betting_history, community_cards_str);
    std::cout << "Subgame root state:\n" << subgame_root.to_string(false) << std::endl;
    
    std::cout << "DEBUG: Subgame root - game_over=" << subgame_root.game_over 
              << ", is_chance=" << subgame_root.is_chance_node()
              << ", current_player=" << subgame_root.current_player()
              << ", round=" << subgame_root.round << std::endl;

    G_ORIGINAL_COMMUNITY_CARDS = subgame_root.community_cards;

    std::vector<Card> known_cards = subgame_root.community_cards;
    std::vector<std::pair<Card, Card>> possible_hands_p0 = get_possible_hole_cards_excluding(known_cards);
    std::vector<std::pair<Card, Card>> possible_hands_p1 = get_possible_hole_cards_excluding(known_cards);
    
    std::cout << "Found " << possible_hands_p0.size() << " possible hands for P0" << std::endl;
    std::cout << "Found " << possible_hands_p1.size() << " possible hands for P1" << std::endl;
    std::cout << "Starting subgame MCCFR for " << num_iterations << " iterations..." << std::endl;

    std::vector<std::thread> threads;
    int iterations_per_thread_base = num_iterations / num_threads_to_use;
    int iterations_remainder = num_iterations % num_threads_to_use;

    auto mccfr_subgame_task = [&](int iterations_for_this_thread, int thread_id) {
        std::mt19937 thread_rng(std::random_device{}() + thread_id);
        
        std::cout << "DEBUG: Thread " << thread_id << " starting with " << iterations_for_this_thread << " iterations" << std::endl;
        
        for (int i = 0; i < iterations_for_this_thread; ++i) {
            if (i < 3) {
                std::cout << "DEBUG: Thread " << thread_id << " iteration " << i << std::endl;
            }
            
            // Sample hole cards for P0
            std::uniform_int_distribution<int> hand_dist_p0(0, possible_hands_p0.size() - 1);
            int sampled_hand_idx_p0 = hand_dist_p0(thread_rng);
            
            if (i < 3) {
                std::cout << "DEBUG: Thread " << thread_id << " sampled P0 hand idx=" << sampled_hand_idx_p0 
                          << " (" << possible_hands_p0[sampled_hand_idx_p0].first.toString() 
                          << " " << possible_hands_p0[sampled_hand_idx_p0].second.toString() << ")" << std::endl;
            }
            
            // Sample hole cards for P1, ensuring no overlap with P0's cards
            std::uniform_int_distribution<int> hand_dist_p1(0, possible_hands_p1.size() - 1);
            int sampled_hand_idx_p1;
            int max_attempts = 100;
            int attempts = 0;
            
            do {
                sampled_hand_idx_p1 = hand_dist_p1(thread_rng);
                attempts++;
            } while (attempts < max_attempts && 
                     hands_overlap(possible_hands_p0[sampled_hand_idx_p0], possible_hands_p1[sampled_hand_idx_p1]));
            
            if (attempts >= max_attempts) {
                std::cerr << "Warning: Could not find non-overlapping hands after " << max_attempts << " attempts" << std::endl;
                continue;
            }
            
            if (i < 3) {
                std::cout << "DEBUG: Thread " << thread_id << " sampled P1 hand idx=" << sampled_hand_idx_p1 
                          << " (" << possible_hands_p1[sampled_hand_idx_p1].first.toString() 
                          << " " << possible_hands_p1[sampled_hand_idx_p1].second.toString() << ")" << std::endl;
            }
            
            // Alternate between updating player 0 and player 1
            int updating_player = (G_CFR_COMPLETED_ITERATIONS_COUNT.load(std::memory_order_relaxed) + i) % 2;
            
            if (i < 3) {
                std::cout << "DEBUG: Thread " << thread_id << " calling mccfr_subgame_recursive for iteration " << i 
                          << " updating player " << updating_player << std::endl;
            }
            
            mccfr_subgame_recursive(subgame_root, updating_player, 
                                   possible_hands_p0, possible_hands_p1,
                                   sampled_hand_idx_p0, sampled_hand_idx_p1, thread_rng);
            
            if (i < 3) {
                std::cout << "DEBUG: Thread " << thread_id << " completed mccfr_subgame_recursive for iteration " << i << std::endl;
            }
            
            G_CFR_COMPLETED_ITERATIONS_COUNT.fetch_add(1, std::memory_order_relaxed);
        }
        
        std::cout << "DEBUG: Thread " << thread_id << " completed all iterations" << std::endl;
    };

    for (int t = 0; t < num_threads_to_use; ++t) {
        int iterations_for_this_thread = iterations_per_thread_base + (t < iterations_remainder ? 1 : 0);
        if (iterations_for_this_thread > 0) {
            threads.emplace_back(mccfr_subgame_task, iterations_for_this_thread, t);
        }
    }

    // Progress monitoring
    int report_interval_ms = 500;
    int iterations_reported_at_last_print = -1;
    bool first_report_triggered = false;

    if (num_iterations > 0) {
        std::cout << "Subgame MCCFR iteration progress (target " << num_iterations << "):" << std::endl;
    }

    while(G_CFR_COMPLETED_ITERATIONS_COUNT.load(std::memory_order_acquire) < num_iterations) {
        std::this_thread::sleep_for(std::chrono::milliseconds(report_interval_ms));
        int current_completed = G_CFR_COMPLETED_ITERATIONS_COUNT.load(std::memory_order_acquire);
        
        if (!first_report_triggered || current_completed > iterations_reported_at_last_print) {
            std::cout << "  Completed " << current_completed << " / " << num_iterations << " iterations." << std::endl;
            iterations_reported_at_last_print = current_completed;
            first_report_triggered = true;
        }
        if (current_completed >= num_iterations) break; 
    }

    // Join threads
    for (auto& th : threads) {
        if (th.joinable()) {
            th.join();
        }
    }
    
    std::cout << "Subgame MCCFR iterations complete. Saving strategy..." << std::endl;
    
    std::cout << "DEBUG: Found " << G_CFR_INFOSET_ACTIONS.size() << " infosets total" << std::endl;
    std::cout << "DEBUG: G_CFR_STRATEGY_SUM size: " << G_CFR_STRATEGY_SUM.size() << std::endl;
    std::cout << "DEBUG: G_CFR_REGRET_SUM size: " << G_CFR_REGRET_SUM.size() << std::endl;
    
    // Save strategy to CSV
    std::ofstream outfile(output_csv_filename);
    outfile << "infoset_key,action_string,probability\n";

    int lines_written = 0;
    for (const auto& pair : G_CFR_INFOSET_ACTIONS) {
        const std::string& infoset_key = pair.first;
        const auto& actions_for_infoset = pair.second;

        if (lines_written < 5) {
            std::cout << "DEBUG: Processing infoset '" << infoset_key << "' with " << actions_for_infoset.size() << " actions" << std::endl;
        }

        if (G_CFR_STRATEGY_SUM.count(infoset_key)) {
            const std::vector<double>& summed_strategy = G_CFR_STRATEGY_SUM.at(infoset_key);
            double total_summed_strategy = 0.0;
            for (double prob_sum : summed_strategy) {
                total_summed_strategy += prob_sum;
            }
            
            if (lines_written < 5) {
                std::cout << "DEBUG: Infoset '" << infoset_key << "' total_summed_strategy=" << total_summed_strategy << std::endl;
            }

            if (total_summed_strategy > 1e-9) {
                for (size_t a = 0; a < summed_strategy.size(); ++a) {
                    double avg_prob = summed_strategy[a] / total_summed_strategy;
                    outfile << "\"" << infoset_key << "\",\"" 
                            << action_to_string(actions_for_infoset[a].first, actions_for_infoset[a].second) 
                            << "\"," << std::fixed << std::setprecision(5) << avg_prob << "\n";
                    lines_written++;
                }
            } else {
                if (!actions_for_infoset.empty()) {
                    double uniform_prob = 1.0 / static_cast<double>(actions_for_infoset.size());
                    for (size_t a = 0; a < actions_for_infoset.size(); ++a) {
                        outfile << "\"" << infoset_key << "\",\""
                                << action_to_string(actions_for_infoset[a].first, actions_for_infoset[a].second)
                                << "\"," << std::fixed << std::setprecision(5) << uniform_prob << "\n";
                        lines_written++;
                    }
                }
            }
        }
    }
    outfile.close();
    
    std::cout << "DEBUG: Wrote " << lines_written << " lines to CSV file" << std::endl;
    std::cout << "Subgame strategy saved to " << output_csv_filename << std::endl;
}

// Main function to solve subgame using MCCFR with exploitability tracking
void solve_subgame_mccfr_with_exploitability_tracking(const HUGame& game_config,
                        const std::string& betting_history,
                        const std::string& community_cards_str,
                        int num_iterations,
                        const std::string& output_csv_filename,
                        int num_threads_to_use = 1,
                        int exploitability_check_interval = 100,
                        int exploitability_simulations = 500) {
    
    // Clear global CFR data structures
    G_CFR_REGRET_SUM.clear();
    G_CFR_STRATEGY_SUM.clear();
    G_CFR_INFOSET_ACTIONS.clear();
    G_CFR_COMPLETED_ITERATIONS_COUNT = 0;

    std::cout << "Initializing subgame state with exploitability tracking..." << std::endl;
    std::cout << "Betting history: " << betting_history << std::endl;
    std::cout << "Community cards: " << community_cards_str << std::endl;
    std::cout << "Exploitability will be computed every " << exploitability_check_interval 
              << " iterations using " << exploitability_simulations << " simulations." << std::endl;

    HUState subgame_root = initialize_subgame_state(game_config, betting_history, community_cards_str);
    std::cout << "Subgame root state:\n" << subgame_root.to_string(false) << std::endl;
    
    std::cout << "DEBUG: Subgame root - game_over=" << subgame_root.game_over 
              << ", is_chance=" << subgame_root.is_chance_node()
              << ", current_player=" << subgame_root.current_player()
              << ", round=" << subgame_root.round << std::endl;

    G_ORIGINAL_COMMUNITY_CARDS = subgame_root.community_cards;

    std::vector<Card> known_cards = subgame_root.community_cards;
    std::vector<std::pair<Card, Card>> possible_hands_p0 = get_possible_hole_cards_excluding(known_cards);
    std::vector<std::pair<Card, Card>> possible_hands_p1 = get_possible_hole_cards_excluding(known_cards);
    
    std::cout << "Found " << possible_hands_p0.size() << " possible hands for P0" << std::endl;
    std::cout << "Found " << possible_hands_p1.size() << " possible hands for P1" << std::endl;
    std::cout << "Starting subgame MCCFR for " << num_iterations << " iterations..." << std::endl;

    std::vector<std::thread> threads;
    int iterations_per_thread_base = num_iterations / num_threads_to_use;
    int iterations_remainder = num_iterations % num_threads_to_use;

    // Exploitability tracking variables
    std::vector<std::tuple<int, double, double>> exploitability_history; // (iteration, exploitability, exploitability_percentage)
    std::mutex exploitability_history_mutex;
    std::atomic<int> next_exploitability_check(exploitability_check_interval);

    auto mccfr_subgame_task = [&](int iterations_for_this_thread, int thread_id) {
        std::mt19937 thread_rng(std::random_device{}() + thread_id);
        
        std::cout << "DEBUG: Thread " << thread_id << " starting with " << iterations_for_this_thread << " iterations" << std::endl;
        
        for (int i = 0; i < iterations_for_this_thread; ++i) {
            if (i < 3) {
                std::cout << "DEBUG: Thread " << thread_id << " iteration " << i << std::endl;
            }
            
            // Sample hole cards for P0
            std::uniform_int_distribution<int> hand_dist_p0(0, possible_hands_p0.size() - 1);
            int sampled_hand_idx_p0 = hand_dist_p0(thread_rng);
            
            // Sample hole cards for P1, ensuring no overlap with P0's cards
            std::uniform_int_distribution<int> hand_dist_p1(0, possible_hands_p1.size() - 1);
            int sampled_hand_idx_p1;
            int max_attempts = 100;
            int attempts = 0;
            
            do {
                sampled_hand_idx_p1 = hand_dist_p1(thread_rng);
                attempts++;
            } while (attempts < max_attempts && 
                     hands_overlap(possible_hands_p0[sampled_hand_idx_p0], possible_hands_p1[sampled_hand_idx_p1]));
            
            if (attempts >= max_attempts) {
                std::cerr << "Warning: Could not find non-overlapping hands after " << max_attempts << " attempts" << std::endl;
                continue;
            }
            
            // Alternate between updating player 0 and player 1
            int updating_player = (G_CFR_COMPLETED_ITERATIONS_COUNT.load(std::memory_order_relaxed) + i) % 2;
            
            mccfr_subgame_recursive(subgame_root, updating_player, 
                                   possible_hands_p0, possible_hands_p1,
                                   sampled_hand_idx_p0, sampled_hand_idx_p1, thread_rng);
            
            int completed_count = G_CFR_COMPLETED_ITERATIONS_COUNT.fetch_add(1, std::memory_order_relaxed) + 1;
            
            // Check if this thread should perform an exploitability check
            int expected_check = next_exploitability_check.load(std::memory_order_acquire);
            if (completed_count >= expected_check && completed_count > 0) {
                // Try to claim this exploitability check
                if (next_exploitability_check.compare_exchange_strong(expected_check, expected_check + exploitability_check_interval)) {
                    std::cout << "\n--- Subgame Exploitability Check at Iteration " << completed_count << " (Thread " << thread_id << ") ---" << std::endl;
                    
                    try {
                        double exploitability = compute_current_strategy_exploitability(subgame_root, exploitability_simulations);
                        
                        {
                            std::lock_guard<std::mutex> lock(exploitability_history_mutex);
                            exploitability_history.push_back(std::make_tuple(completed_count, exploitability, 0.0));
                        }
                        
                        std::cout << "Iteration " << completed_count << " - Subgame Exploitability: " << std::fixed 
                                  << std::setprecision(6) << exploitability << std::endl;
                        
                        // Show improvement if we have previous data
                        {
                            std::lock_guard<std::mutex> lock(exploitability_history_mutex);
                            if (exploitability_history.size() > 1) {
                                double prev_exploitability = std::get<1>(exploitability_history[exploitability_history.size() - 2]);
                                double improvement = prev_exploitability - exploitability;
                                std::cout << "Improvement since last check: " << std::fixed 
                                          << std::setprecision(6) << improvement;
                                if (prev_exploitability > 0) {
                                    double improvement_percent = (improvement / prev_exploitability) * 100.0;
                                    std::cout << " (" << std::fixed << std::setprecision(2) << improvement_percent << "% reduction)";
                                }
                                std::cout << std::endl;
                            }
                        }
                        std::cout << "--- End Subgame Exploitability Check ---\n" << std::endl;
                        
                    } catch (const std::exception& e) {
                        std::cerr << "Error computing subgame exploitability: " << e.what() << std::endl;
                    }
                }
            }
        }
        
        std::cout << "DEBUG: Thread " << thread_id << " completed all iterations" << std::endl;
    };

    for (int t = 0; t < num_threads_to_use; ++t) {
        int iterations_for_this_thread = iterations_per_thread_base + (t < iterations_remainder ? 1 : 0);
        if (iterations_for_this_thread > 0) {
            threads.emplace_back(mccfr_subgame_task, iterations_for_this_thread, t);
        }
    }

    // Progress monitoring with exploitability tracking
    int report_interval_ms = 500;
    int iterations_reported_at_last_print = -1;
    bool first_report_triggered = false;
    int last_exploitability_check_iteration = -1;

    if (num_iterations > 0) {
        std::cout << "Subgame MCCFR iteration progress with exploitability tracking (target " << num_iterations << "):" << std::endl;
    }

    while(G_CFR_COMPLETED_ITERATIONS_COUNT.load(std::memory_order_acquire) < num_iterations) {
        std::this_thread::sleep_for(std::chrono::milliseconds(report_interval_ms));
        int current_completed = G_CFR_COMPLETED_ITERATIONS_COUNT.load(std::memory_order_acquire);
        
        // Check if we should compute exploitability
        if (current_completed > 0 && 
            current_completed >= last_exploitability_check_iteration + exploitability_check_interval &&
            current_completed != last_exploitability_check_iteration) {
            
            last_exploitability_check_iteration = current_completed;
            
            std::cout << "\n--- Subgame Exploitability Check at Iteration " << current_completed << " ---" << std::endl;
            
            try {
                double exploitability = compute_current_strategy_exploitability(subgame_root, exploitability_simulations);
                exploitability_history.push_back(std::make_tuple(current_completed, exploitability, 0.0));
                
                std::cout << "Iteration " << current_completed << " - Subgame Exploitability: " << std::fixed 
                          << std::setprecision(6) << exploitability << std::endl;
                
                // Show improvement if we have previous data
                if (exploitability_history.size() > 1) {
                    double prev_exploitability = std::get<1>(exploitability_history[exploitability_history.size() - 2]);
                    double improvement = prev_exploitability - exploitability;
                    std::cout << "Improvement since last check: " << std::fixed 
                              << std::setprecision(6) << improvement;
                    if (prev_exploitability > 0) {
                        double improvement_percent = (improvement / prev_exploitability) * 100.0;
                        std::cout << " (" << std::fixed << std::setprecision(2) << improvement_percent << "% reduction)";
                    }
                    std::cout << std::endl;
                }
                std::cout << "--- End Subgame Exploitability Check ---\n" << std::endl;
                
            } catch (const std::exception& e) {
                std::cerr << "Error computing subgame exploitability: " << e.what() << std::endl;
            }
        }
        
        if (!first_report_triggered || current_completed > iterations_reported_at_last_print) {
            std::cout << "  Completed " << current_completed << " / " << num_iterations << " iterations." << std::endl;
            iterations_reported_at_last_print = current_completed;
            first_report_triggered = true;
        }
        if (current_completed >= num_iterations) break; 
    }

    // Join threads
    for (auto& th : threads) {
        if (th.joinable()) {
            th.join();
        }
    }
    
    // Final completion message
    int final_completed_count = G_CFR_COMPLETED_ITERATIONS_COUNT.load(std::memory_order_acquire);
    if (num_iterations > 0) {
        if (final_completed_count > iterations_reported_at_last_print || !first_report_triggered) {
             std::cout << "  Completed " << final_completed_count << " / " << num_iterations << " iterations." << std::endl;
        }
        std::cout << "Subgame MCCFR Iteration " << final_completed_count << "/" << num_iterations << " (All threads complete)" << std::endl;
    }

    // Final exploitability computation
    std::cout << "\nComputing final subgame exploitability..." << std::endl;
    try {
        double final_exploitability = compute_current_strategy_exploitability(subgame_root, exploitability_simulations * 2);
        {
            std::lock_guard<std::mutex> lock(exploitability_history_mutex);
            exploitability_history.push_back(std::make_tuple(final_completed_count, final_exploitability, 0.0));
        }
        std::cout << "Final subgame exploitability: " << std::fixed << std::setprecision(6) << final_exploitability << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Error computing final subgame exploitability: " << e.what() << std::endl;
    }

    // Print exploitability history summary
    {
        std::lock_guard<std::mutex> lock(exploitability_history_mutex);
        if (!exploitability_history.empty()) {
            std::cout << "\nSubgame Exploitability Training History:" << std::endl;
            std::cout << "Iteration\tExploitability\tImprovement" << std::endl;
            for (size_t i = 0; i < exploitability_history.size(); ++i) {
                const auto& entry = exploitability_history[i];
                std::cout << std::get<0>(entry) << "\t\t" << std::fixed << std::setprecision(6) << std::get<1>(entry);
                if (i > 0) {
                    double improvement = std::get<1>(exploitability_history[i-1]) - std::get<1>(entry);
                    std::cout << "\t\t" << std::fixed << std::setprecision(6) << improvement;
                } else {
                    std::cout << "\t\t--";
                }
                std::cout << std::endl;
            }
            
            if (exploitability_history.size() > 1) {
                double initial_exp = std::get<1>(exploitability_history[0]);
                double final_exp = std::get<1>(exploitability_history.back());
                double total_improvement = initial_exp - final_exp;
                std::cout << "Total exploitability improvement: " << std::fixed << std::setprecision(6) 
                          << total_improvement << " (" << std::fixed << std::setprecision(2) 
                          << (total_improvement / initial_exp * 100.0) << "% reduction)" << std::endl;
            }
        }
    }
    
    std::cout << "Subgame MCCFR iterations complete. Saving strategy..." << std::endl;
    
    std::cout << "DEBUG: Found " << G_CFR_INFOSET_ACTIONS.size() << " infosets total" << std::endl;
    std::cout << "DEBUG: G_CFR_STRATEGY_SUM size: " << G_CFR_STRATEGY_SUM.size() << std::endl;
    std::cout << "DEBUG: G_CFR_REGRET_SUM size: " << G_CFR_REGRET_SUM.size() << std::endl;
    
    // Save strategy to CSV
    std::ofstream outfile(output_csv_filename);
    outfile << "infoset_key,action_string,probability\n";

    int lines_written = 0;
    for (const auto& pair : G_CFR_INFOSET_ACTIONS) {
        const std::string& infoset_key = pair.first;
        const auto& actions_for_infoset = pair.second;

        if (G_CFR_STRATEGY_SUM.count(infoset_key)) {
            const std::vector<double>& summed_strategy = G_CFR_STRATEGY_SUM.at(infoset_key);
            double total_summed_strategy = 0.0;
            for (double prob_sum : summed_strategy) {
                total_summed_strategy += prob_sum;
            }

            if (total_summed_strategy > 1e-9) {
                for (size_t a = 0; a < summed_strategy.size(); ++a) {
                    double avg_prob = summed_strategy[a] / total_summed_strategy;
                    outfile << "\"" << infoset_key << "\",\"" 
                            << action_to_string(actions_for_infoset[a].first, actions_for_infoset[a].second) 
                            << "\"," << std::fixed << std::setprecision(5) << avg_prob << "\n";
                    lines_written++;
                }
            } else {
                if (!actions_for_infoset.empty()) {
                    double uniform_prob = 1.0 / static_cast<double>(actions_for_infoset.size());
                    for (size_t a = 0; a < actions_for_infoset.size(); ++a) {
                        outfile << "\"" << infoset_key << "\",\""
                                << action_to_string(actions_for_infoset[a].first, actions_for_infoset[a].second)
                                << "\"," << std::fixed << std::setprecision(5) << uniform_prob << "\n";
                        lines_written++;
                    }
                }
            }
        }
    }
    outfile.close();
    
    // Save exploitability history to file
    std::string exploitability_log_filename = output_csv_filename + "_subgame_exploitability_log.csv";
    std::ofstream exp_outfile(exploitability_log_filename);
    exp_outfile << "iteration,exploitability\n";
    for (const auto& entry : exploitability_history) {
        exp_outfile << std::get<0>(entry) << "," << std::fixed << std::setprecision(6) << std::get<1>(entry) << "\n";
    }
    exp_outfile.close();
    
    std::cout << "DEBUG: Wrote " << lines_written << " lines to CSV file" << std::endl;
    std::cout << "Subgame strategy saved to " << output_csv_filename << std::endl;
    std::cout << "Subgame exploitability history saved to " << exploitability_log_filename << std::endl;
}

int main(int argc, char* argv[]) {
    // Attempt to preload clusters at the start of the program.
    // Ensure the CSV files are in the correct path (e.g., ./utils/) relative to where the executable is run.
    preloadClusters_ic();

    // Default max actions per simulation to prevent infinite loops in manual or auto play.
    int max_actions_per_sim_arg = 1e9; 
    // initial_num_simulations_arg is no longer used for a global simulation,
    // but parsing logic is kept for max_actions_per_sim_arg.
    int initial_num_simulations_arg = 1e9; 


    if (argc > 1) {
        try {
            // The first argument might be intended for max_actions if only one is provided,
            // or num_simulations if two are provided.
            // For simplicity, let's assume argv[1] is num_simulations (now unused)
            // and argv[2] is max_actions.
            // If you want to control these via args, adjust parsing.
            // For now, we'll primarily use max_actions_per_sim_arg.
            initial_num_simulations_arg = std::stoi(argv[1]); // Not directly used in interactive mode
        } catch (const std::exception& e) {
            std::cerr << "Warning: Invalid argument for initial number of simulations: " << argv[1]
                      << ". This argument is not directly used in interactive mode." << std::endl;
        }
    }
    if (argc > 2) {
        try {
            max_actions_per_sim_arg = std::stoi(argv[2]);
        } catch (const std::exception& e) {
            std::cerr << "Invalid argument for max actions per simulation: " << argv[2]
                      << ". Using default: " << max_actions_per_sim_arg << std::endl;
        }
    }

    std::cout << "=== HEADS-UP POKER MCCFR SOLVER ===" << std::endl;
    std::cout << "This solver allows you to:" << std::endl;
    std::cout << "- Configure game parameters" << std::endl;
    std::cout << "- Play poker interactively" << std::endl;
    std::cout << "- Run MCCFR solving from any state" << std::endl;
    std::cout << "- Solve subgames using specialized algorithms" << std::endl;
    std::cout << "- Count information sets via Monte Carlo" << std::endl;
    std::cout << std::endl;

    // --- User-configurable game settings ---
    double user_initial_stack = 100.0;
    double user_small_blind = 0.5;
    double user_big_blind = 1.0;
    float user_rake_percentage = 0.07;
    std::vector<double> user_pot_fraction_bet_sizes = {0.33, 0.5, 0.75, 1.0};
    std::vector<double> user_fixed_bet_sizes_bb = {2.5, 3.0, 4.0};
    double user_all_in_threshold = 1.5;        // Default 150% (1.5)
    double user_force_all_in_threshold = 0.2;  // Default 20% (0.2)
    double user_merging_threshold = 0.1;       // Default 10% (0.1)

    std::string line_input;

    std::cout << "--- Configure Game Parameters ---" << std::endl;

    // Initial Stack
    std::cout << "Enter initial stack for each player (e.g., 100.0): ";
    while (!(std::cin >> user_initial_stack) || user_initial_stack <= 0) {
        std::cout << "Invalid input. Please enter a positive number for initial stack: ";
        std::cin.clear();
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    // Small Blind
    std::cout << "Enter small blind amount (e.g., 0.5): ";
    while (!(std::cin >> user_small_blind) || user_small_blind <= 0) {
        std::cout << "Invalid input. Please enter a positive number for small blind: ";
        std::cin.clear();
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    // Big Blind
    std::cout << "Enter big blind amount (e.g., 1.0 - must be >= small blind): ";
    while (!(std::cin >> user_big_blind) || user_big_blind < user_small_blind) {
        std::cout << "Invalid input. Please enter a positive number for big blind (>= small blind): ";
        std::cin.clear();
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    // Rake Percentage
    std::cout << "Enter rake percentage (e.g., 0.07 for 7%): ";
    while (!(std::cin >> user_rake_percentage) || user_rake_percentage < 0 || user_rake_percentage >= 1) {
        std::cout << "Invalid input. Please enter a number between 0.0 and 0.99 for rake percentage: ";
        std::cin.clear();
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    // All-in Threshold
    std::cout << "Enter all-in threshold (e.g., 1.5 for 150% of pot): ";
    while (!(std::cin >> user_all_in_threshold) || user_all_in_threshold <= 0) {
        std::cout << "Invalid input. Please enter a positive number for all-in threshold: ";
        std::cin.clear();
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    // Force All-in Threshold
    std::cout << "Enter force all-in threshold (e.g., 0.2 for 20% of stack): ";
    while (!(std::cin >> user_force_all_in_threshold) || user_force_all_in_threshold < 0 || user_force_all_in_threshold > 1) {
        std::cout << "Invalid input. Please enter a number between 0.0 and 1.0 for force all-in threshold: ";
        std::cin.clear();
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    // Merging Threshold
    std::cout << "Enter merging threshold (e.g., 0.1 for 10%): ";
    while (!(std::cin >> user_merging_threshold) || user_merging_threshold < 0 || user_merging_threshold > 1) {
        std::cout << "Invalid input. Please enter a number between 0.0 and 1.0 for merging threshold: ";
        std::cin.clear();
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    // Pot Fraction Bet Sizes
    std::cout << "Enter space-separated pot fraction bet sizes (e.g., 0.33 0.5 0.75 1.0): ";
    std::getline(std::cin, line_input);
    std::istringstream iss_pot_frac(line_input);
    user_pot_fraction_bet_sizes.clear();
    double val;
    while (iss_pot_frac >> val) {
        if (val > 0) user_pot_fraction_bet_sizes.push_back(val);
    }
    if (user_pot_fraction_bet_sizes.empty()) {
        std::cout << "No valid pot fraction bet sizes entered. Using default: 0.33 0.5 0.75 1.0" << std::endl;
        user_pot_fraction_bet_sizes = {0.33, 0.5, 0.75, 1.0};
    }

    // Fixed BB Bet Sizes
    std::cout << "Enter space-separated fixed big blind bet sizes (e.g., 2.5 3.0 4.0): ";
    std::getline(std::cin, line_input);
    std::istringstream iss_fixed_bb(line_input);
    user_fixed_bet_sizes_bb.clear();
    while (iss_fixed_bb >> val) {
        if (val > 0) user_fixed_bet_sizes_bb.push_back(val);
    }
    if (user_fixed_bet_sizes_bb.empty()) {
        std::cout << "No valid fixed BB bet sizes entered. Using default: 2.5 3.0 4.0" << std::endl;
        user_fixed_bet_sizes_bb = {2.5, 3.0, 4.0};
    }
    
    std::cout << "--- Game configuration complete ---" << std::endl;

    HUGame game_config(
        2, // num_players fixed at 2 for Heads-Up
        user_initial_stack,
        user_small_blind,
        user_big_blind,
        user_rake_percentage,
        user_pot_fraction_bet_sizes,
        user_fixed_bet_sizes_bb,
        user_all_in_threshold,
        user_force_all_in_threshold,
        user_merging_threshold
    );

    std::cout << "\nStarting Interactive Heads-Up Poker Game..." << std::endl;
    std::cout << "Game settings: Initial Stack=" << game_config.initial_stack
              << ", SB=" << game_config.small_blind_amount
              << ", BB=" << game_config.big_blind_amount
              << ", Rake=" << (game_config.rake_percentage * 100) << "%" << std::endl;
    std::cout << "Thresholds: All-in=" << user_all_in_threshold
              << ", Force All-in=" << user_force_all_in_threshold
              << ", Merging=" << user_merging_threshold << std::endl;
    std::cout << "Pot fraction bet sizes: ";
    for (size_t i = 0; i < game_config.pot_fraction_bet_sizes.size(); ++i) {
        if (i > 0) std::cout << ", ";
        std::cout << game_config.pot_fraction_bet_sizes[i];
    }
    std::cout << std::endl;
    std::cout << "Fixed BB bet sizes: ";
    for (size_t i = 0; i < game_config.fixed_bet_sizes_bb.size(); ++i) {
        if (i > 0) std::cout << ", ";
        std::cout << game_config.fixed_bet_sizes_bb[i];
    }
    std::cout << std::endl;
    std::cout << "Max actions for manual play / per simulation: " << max_actions_per_sim_arg << std::endl;

    HUState current_game_state = game_config.new_initial_state();
    int actions_taken_manual = 0;

    while (!current_game_state.game_over && actions_taken_manual < max_actions_per_sim_arg) {
        actions_taken_manual++;
        std::cout << "\n-----------------------------------\n";
        std::cout << "MANUAL PLAY - Action #" << actions_taken_manual << std::endl;
        // Display full state including cards for the current player for decision making
        std::cout << current_game_state.to_string(true) << std::endl;

        if (current_game_state.is_chance_node()) {
            std::cout << "Dealing cards (automatic action)...\n";
            current_game_state.apply_action(Action::DEAL);
            // Show state after dealing
            std::cout << "State after dealing:\n" << current_game_state.to_string(true) << std::endl;
            std::cout << "Press Enter to continue...";
            std::string dummy_line;
            std::getline(std::cin, dummy_line); // Single getline is sufficient
            continue;
        }

        std::vector<std::pair<Action, double>> legal_actions = current_game_state.legal_actions(
            current_game_state.game.pot_fraction_bet_sizes,
            current_game_state.game.fixed_bet_sizes_bb
        );

        if (legal_actions.empty() && !current_game_state.game_over) {
            std::cout << "No legal actions available for player " << current_game_state.current_player()
                      << ", but game not over. This might indicate an all-in situation or end of round." << std::endl;
            std::cout << "The game might auto-progress or end soon." << std::endl;
            // This state might occur if a player is all-in and it's their turn technically,
            // but they have no more actions. The game should proceed.
            // For now, we break the manual loop; the game state will be printed.
            // Ideally, the game logic should advance past such states automatically if no user input is required.
            // However, apply_action is what advances state, so if legal_actions is empty, we can't call it.
            // This implies the game should have already transitioned or ended.
            break;
        }
        
        std::cout << "\nPlayer " << current_game_state.current_player() << "'s turn. Hole Cards: "
                  << current_game_state.get_player_hole_cards_str(current_game_state.current_player(), true)
                  << ". Legal actions:\n";

        for (size_t i = 0; i < legal_actions.size(); ++i) {
            std::cout << i << ": " << action_to_string(legal_actions[i].first, legal_actions[i].second) << "\n";
        }
        std::cout << legal_actions.size() << ": Simulate infosets from this state\n";
        std::cout << legal_actions.size() + 1 << ": Solve from this state (MCCFR)\n";
        std::cout << legal_actions.size() + 2 << ": Solve subgame with exploitability tracking (Subgame MCCFR)\n";
        std::cout << legal_actions.size() + 3 << ": Compute exploitability of current CFR strategy (requires previous MCCFR run)\n";
        std::cout << legal_actions.size() + 4 << ": Compute exploitability from strategy CSV file (MES formula)\n";
        std::cout << legal_actions.size() + 5 << ": Solve with exploitability tracking (MCCFR)\n";
        std::cout << legal_actions.size() + 6 << ": Solve from this state (ACCELERATED MCCFR) - NEW!\n";
        std::cout << legal_actions.size() + 7 << ": Quit game\n";

        int choice_val;
        std::cout << "Enter the number of your chosen action: ";
        while (!(std::cin >> choice_val) || choice_val < 0 || choice_val > static_cast<int>(legal_actions.size() + 7)) {
            std::cout << "Invalid input. Please enter a number between 0 and " << legal_actions.size() + 7 << ": ";
            std::cin.clear();
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        }
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // Consume the newline character

        if (choice_val == static_cast<int>(legal_actions.size())) {
            // Simulate infosets
            int num_sims_for_subtree;
            std::cout << "Enter number of simulations to run from this state: ";
            while (!(std::cin >> num_sims_for_subtree) || num_sims_for_subtree <= 0) {
                std::cout << "Invalid input. Please enter a positive integer: ";
                std::cin.clear();
                std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            }
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // Consume newline

            int infoset_count_subtree = count_unique_infosets_mc(
                current_game_state, // Pass current game state
                num_sims_for_subtree,
                max_actions_per_sim_arg
            );
            std::cout << "\nEstimated number of unique player information sets from current state: "
                      << infoset_count_subtree << std::endl;
            std::cout << "Continuing manual play from the same state...\n";
            actions_taken_manual--; // Don't count simulation as a game action for the limit
        } else if (choice_val == static_cast<int>(legal_actions.size() + 1)) {
            // Solve from this state (MCCFR)
            int num_mccfr_iterations;
            std::cout << "Enter number of MCCFR iterations to run from this state: ";
            while (!(std::cin >> num_mccfr_iterations) || num_mccfr_iterations <= 0) {
                std::cout << "Invalid input. Please enter a positive integer for iterations: ";
                std::cin.clear();
                std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            }
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // Consume newline
            
            unsigned int max_hw_threads = std::thread::hardware_concurrency();
            if (max_hw_threads == 0) {
                std::cout << "Warning: Could not determine maximum hardware threads. Defaulting to 1." << std::endl;
                max_hw_threads = 1;
            }
            std::cout << "Your system suggests up to " << max_hw_threads << " concurrent threads." << std::endl;
            
            int num_threads_to_use = 1;
            if (max_hw_threads > 1) {
                std::cout << "Enter number of threads to use for MCCFR (1 to " << max_hw_threads << "): ";
                while (!(std::cin >> num_threads_to_use) || num_threads_to_use < 1 || num_threads_to_use > static_cast<int>(max_hw_threads)) {
                    std::cout << "Invalid input. Please enter a number between 1 and " << max_hw_threads << ": ";
                    std::cin.clear();
                    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                }
                std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // Consume newline
            } else {
                std::cout << "Using 1 thread for MCCFR." << std::endl;
                num_threads_to_use = 1;
            }
            
            std::string mccfr_csv_filename = "mccfr_strategy_from_state.csv";
            solve_from_state_mccfr(current_game_state, num_mccfr_iterations, mccfr_csv_filename, num_threads_to_use);
            
            std::cout << "MCCFR solving complete. Strategy saved. Continuing manual play from the same state...\n";
            actions_taken_manual--; // Don't count MCCFR solving as a game action for the limit

        } else if (choice_val == static_cast<int>(legal_actions.size() + 2)) {
            // Solve subgame from current state (Subgame MCCFR with Exploitability Tracking)
            
            // Automatically extract betting history and community cards from current state
            std::string betting_history_input = extract_betting_history_from_state(current_game_state);
            std::string community_cards_input;
            std::cout << "Enter community cards (e.g., 'As Kd Qh' for flop, or empty for preflop): ";
            std::getline(std::cin, community_cards_input);
            
            // Remove trailing whitespace from community cards string
            community_cards_input.erase(community_cards_input.find_last_not_of(" \n\r\t") + 1);
            
            std::cout << "Automatically extracted from current state:" << std::endl;
            std::cout << "  Betting history: '" << betting_history_input << "'" << std::endl;
            std::cout << "  Community cards: '" << community_cards_input << "'" << std::endl;
            std::cout << "  Current round: " << current_game_state.round << std::endl;
            
            int num_subgame_iterations;
            std::cout << "Enter number of Subgame MCCFR iterations: ";
            while (!(std::cin >> num_subgame_iterations) || num_subgame_iterations <= 0) {
                std::cout << "Invalid input. Please enter a positive integer for iterations: ";
                std::cin.clear();
                std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            }
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // Consume newline
            
            // Ask for exploitability check interval (leap parameter)
            int exploitability_check_interval;
            std::cout << "Enter exploitability check interval (every X iterations, e.g., 10 or 100): ";
            while (!(std::cin >> exploitability_check_interval) || exploitability_check_interval <= 0) {
                std::cout << "Invalid input. Please enter a positive integer for check interval: ";
                std::cin.clear();
                std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            }
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // Consume newline
            
            // Ask for number of simulations for exploitability calculation
            int exploitability_simulations;
            std::cout << "Enter number of simulations for each exploitability calculation (e.g., 500): ";
            while (!(std::cin >> exploitability_simulations) || exploitability_simulations <= 0) {
                std::cout << "Invalid input. Please enter a positive integer for simulations: ";
                std::cin.clear();
                std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            }
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // Consume newline
            
            unsigned int max_hw_threads = std::thread::hardware_concurrency();
            if (max_hw_threads == 0) {
                std::cout << "Warning: Could not determine maximum hardware threads. Defaulting to 1." << std::endl;
                max_hw_threads = 1;
            }
            std::cout << "Your system suggests up to " << max_hw_threads << " concurrent threads." << std::endl;
            
            int num_threads_to_use = 1;
            if (max_hw_threads > 1) {
                std::cout << "Enter number of threads to use for Subgame MCCFR (1 to " << max_hw_threads << "): ";
                while (!(std::cin >> num_threads_to_use) || num_threads_to_use < 1 || num_threads_to_use > static_cast<int>(max_hw_threads)) {
                    std::cout << "Invalid input. Please enter a number between 1 and " << max_hw_threads << ": ";
                    std::cin.clear();
                    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                }
                std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // Consume newline
            } else {
                std::cout << "Using 1 thread for Subgame MCCFR." << std::endl;
                num_threads_to_use = 1;
            }
            
            std::string subgame_mccfr_csv_filename = "subgame_mccfr_strategy_from_state.csv";
            solve_subgame_mccfr_with_exploitability_tracking(game_config, betting_history_input, community_cards_input, 
                                                           num_subgame_iterations, subgame_mccfr_csv_filename, 
                                                           num_threads_to_use, exploitability_check_interval, 
                                                           exploitability_simulations);
            
            std::cout << "Subgame MCCFR solving with exploitability tracking complete. Strategy saved. Continuing manual play from the same state...\n";
            actions_taken_manual--; // Don't count Subgame MCCFR solving as a game action for the limit

        } else if (choice_val == static_cast<int>(legal_actions.size() + 3)) {
            // Compute exploitability of current CFR strategy
            std::cout << "\n=== Computing Exploitability of Current CFR Strategy ===" << std::endl;
            std::cout << "This uses the MES (Maximally Exploitative Strategy) formula:" << std::endl;
            std::cout << "- For Non-raked games: (MES_EV[P0] + MES_EV[P1]) × 0.5" << std::endl;
            std::cout << "- For Raked games: ((MES_EV[P0] - Current_EV[P0]) + (MES_EV[P1] - Current_EV[P1])) × 0.5" << std::endl;
            std::cout << "Note: This requires that MCCFR has been run previously to generate strategy data." << std::endl;
            
            int num_simulations;
            std::cout << "Enter number of simulations for exploitability calculation: ";
            while (!(std::cin >> num_simulations) || num_simulations <= 0) {
                std::cout << "Invalid input. Please enter a positive integer: ";
                std::cin.clear();
                std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            }
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // Consume newline
            
            double exploitability = compute_current_strategy_exploitability(current_game_state, num_simulations);
            std::cout << "\n=== Exploitability Result ===" << std::endl;
            std::cout << "Final Exploitability: " << std::fixed << std::setprecision(6) << exploitability << std::endl;
            std::cout << "=========================" << std::endl;
            actions_taken_manual--; // Don't count exploitability calculation as a game action for the limit

        } else if (choice_val == static_cast<int>(legal_actions.size() + 4)) {
            // Compute exploitability from strategy file
            std::cout << "\n=== Computing Exploitability from Strategy CSV File ===" << std::endl;
            std::cout << "This loads a strategy from a CSV file and computes its exploitability using the MES formula:" << std::endl;
            std::cout << "- For Non-raked games: (MES_EV[P0] + MES_EV[P1]) × 0.5" << std::endl;
            std::cout << "- For Raked games: ((MES_EV[P0] - Current_EV[P0]) + (MES_EV[P1] - Current_EV[P1])) × 0.5" << std::endl;
            std::cout << "CSV format expected: infoset_key,action_string,probability" << std::endl;
            
            std::string strategy_csv_filename;
            std::cout << "Enter path to strategy CSV file: ";
            std::getline(std::cin, strategy_csv_filename);
            
            int num_simulations;
            std::cout << "Enter number of simulations for exploitability calculation: ";
            while (!(std::cin >> num_simulations) || num_simulations <= 0) {
                std::cout << "Invalid input. Please enter a positive integer: ";
                std::cin.clear();
                std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            }
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // Consume newline
            
            double exploitability = compute_strategy_file_exploitability(current_game_state, strategy_csv_filename, num_simulations);
            std::cout << "\n=== Exploitability Result ===" << std::endl;
            std::cout << "Final Exploitability: " << std::fixed << std::setprecision(6) << exploitability << std::endl;
            std::cout << "=========================" << std::endl;
            actions_taken_manual--; // Don't count exploitability calculation as a game action for the limit

        } else if (choice_val == static_cast<int>(legal_actions.size() + 5)) {
            // Solve with exploitability tracking (MCCFR)
            int num_mccfr_iterations;
            std::cout << "Enter number of MCCFR iterations to run with exploitability tracking: ";
            while (!(std::cin >> num_mccfr_iterations) || num_mccfr_iterations <= 0) {
                std::cout << "Invalid input. Please enter a positive integer for iterations: ";
                std::cin.clear();
                std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            }
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // Consume newline
            
            unsigned int max_hw_threads = std::thread::hardware_concurrency();
            if (max_hw_threads == 0) {
                std::cout << "Warning: Could not determine maximum hardware threads. Defaulting to 1." << std::endl;
                max_hw_threads = 1;
            }
            std::cout << "Your system suggests up to " << max_hw_threads << " concurrent threads." << std::endl;
            
            int num_threads_to_use = 1;
            if (max_hw_threads > 1) {
                std::cout << "Enter number of threads to use for MCCFR (1 to " << max_hw_threads << "): ";
                while (!(std::cin >> num_threads_to_use) || num_threads_to_use < 1 || num_threads_to_use > static_cast<int>(max_hw_threads)) {
                    std::cout << "Invalid input. Please enter a number between 1 and " << max_hw_threads << ": ";
                    std::cin.clear();
                    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                }
                std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // Consume newline
            } else {
                std::cout << "Using 1 thread for MCCFR." << std::endl;
                num_threads_to_use = 1;
            }
            
            std::string mccfr_csv_filename = "mccfr_strategy_from_state.csv";
            solve_from_state_mccfr(current_game_state, num_mccfr_iterations, mccfr_csv_filename, num_threads_to_use);
            
            std::cout << "MCCFR solving complete. Strategy saved. Continuing manual play from the same state...\n";
            actions_taken_manual--; // Don't count MCCFR solving as a game action for the limit

        } else if (choice_val == static_cast<int>(legal_actions.size() + 6)) {
            // Solve from this state (ACCELERATED MCCFR)
            int num_accelerated_iterations;
            std::cout << "Enter number of ACCELERATED MCCFR iterations: ";
            while (!(std::cin >> num_accelerated_iterations) || num_accelerated_iterations <= 0) {
                std::cout << "Invalid input. Please enter a positive integer for iterations: ";
                std::cin.clear();
                std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            }
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // Consume newline
            
            solve_from_state_mccfr_accelerated(current_game_state, num_accelerated_iterations, "accelerated_mccfr_strategy.csv", 1, false, true, 0);
            
            std::cout << "ACCELERATED MCCFR solving complete. Strategy saved. Continuing manual play from the same state...\n";
            actions_taken_manual--; // Don't count ACCELERATED MCCFR solving as a game action for the limit

        } else if (choice_val == static_cast<int>(legal_actions.size() + 7)) {
            std::cout << "Quitting manual game." << std::endl;
            break;
        } else {
            // Apply chosen game action
            Action chosen_action = legal_actions[choice_val].first;
            double chosen_amount = legal_actions[choice_val].second;

            std::cout << "Applying action: " << action_to_string(chosen_action, chosen_amount) << std::endl;
            try {
                current_game_state.apply_action(chosen_action, chosen_amount);
            } catch (const std::runtime_error& e) {
                std::cerr << "Error applying action: " << e.what() << std::endl;
                std::cerr << "Current state before error:\n" << current_game_state.to_string(true) << std::endl;
                break; // Exit on error
            }
        }
    }

    std::cout << "\n===================================\n";
    if (current_game_state.game_over) {
        std::cout << "Game Over!\n";
    } else if (actions_taken_manual >= max_actions_per_sim_arg) {
        std::cout << "Max actions for manual play reached (" << max_actions_per_sim_arg << ").\n";
    }
    std::cout << "Final State:\n" << current_game_state.to_string(true) << std::endl;

    if(current_game_state.game_over){
        std::vector<double> final_returns = current_game_state.returns();
        std::cout << "Returns: ";
        bool first_return = true;
        for(size_t p_idx = 0; p_idx < final_returns.size(); ++p_idx) {
            if (!first_return) std::cout << ", ";
            std::cout << "P" << p_idx << ": " << final_returns[p_idx];
            first_return = false;
        }
        std::cout << std::endl;
    }
    std::cout << "Exiting infoset counter." << std::endl;

    return 0;
}

// ===============================
// EXPLOITABILITY CALCULATION CODE
// ===============================

// Data structure to store a complete strategy profile
struct StrategyProfile {
    std::map<std::string, std::vector<double>> strategies; // infoset_key -> action probabilities
    std::map<std::string, std::vector<std::pair<Action, double>>> infoset_actions; // infoset_key -> legal actions
};

// Load strategy from CSV file
StrategyProfile load_strategy_from_csv(const std::string& csv_filename) {
    StrategyProfile profile;
    std::ifstream file(csv_filename);
    
    if (!file.is_open()) {
        std::cerr << "Error: Could not open strategy file " << csv_filename << std::endl;
        return profile;
    }
    
    std::string line;
    std::getline(file, line); // Skip header
    
    while (std::getline(file, line)) {
        // Parse CSV line: "infoset_key","action_string",probability
        std::vector<std::string> fields;
        std::string field;
        bool inQuotes = false;
        
        for (char c : line) {
            if (c == '"') {
                inQuotes = !inQuotes;
            } else if (c == ',' && !inQuotes) {
                fields.push_back(field);
                field.clear();
            } else if (c != '"') { // Skip quotes
                field += c;
            }
        }
        fields.push_back(field);
        
        if (fields.size() >= 3) {
            std::string infoset_key = fields[0];
            std::string action_string = fields[1];
            double probability = std::stod(fields[2]);
            
            // Parse action from action_string
            Action action;
            double amount = 0.0;
            
            if (action_string == "CHECK") {
                action = Action::CHECK;
            } else if (action_string == "CALL") {
                action = Action::CALL;
            } else if (action_string == "FOLD") {
                action = Action::FOLD;
            } else if (action_string.find("BET_RAISE") != std::string::npos) {
                action = Action::BET_RAISE;
                // Extract amount from string if needed
                size_t pos = action_string.find("(");
                if (pos != std::string::npos) {
                    size_t end_pos = action_string.find(")", pos);
                    if (end_pos != std::string::npos) {
                        amount = std::stod(action_string.substr(pos + 1, end_pos - pos - 1));
                    }
                }
            } else if (action_string == "ALL_IN") {
                action = Action::ALL_IN;
            }
            
            // Add to strategy profile
            if (profile.strategies.find(infoset_key) == profile.strategies.end()) {
                profile.strategies[infoset_key] = std::vector<double>();
                profile.infoset_actions[infoset_key] = std::vector<std::pair<Action, double>>();
            }
            
            profile.strategies[infoset_key].push_back(probability);
            profile.infoset_actions[infoset_key].push_back({action, amount});
        }
    }
    
    file.close();
    std::cout << "Loaded strategy with " << profile.strategies.size() << " infosets from " << csv_filename << std::endl;
    return profile;
}

// Convert current CFR data to StrategyProfile
StrategyProfile get_current_strategy_profile() {
    StrategyProfile profile;
    
    std::lock_guard<std::mutex> actions_lock(G_CFR_INFOSET_ACTIONS_MUTEX);
    std::lock_guard<std::mutex> strategy_lock(G_CFR_STRATEGY_SUM_MUTEX);
    
    for (const auto& pair : G_CFR_INFOSET_ACTIONS) {
        const std::string& infoset_key = pair.first;
        const auto& actions_for_infoset = pair.second;
        
        profile.infoset_actions[infoset_key] = actions_for_infoset;
        
        if (G_CFR_STRATEGY_SUM.count(infoset_key)) {
            const std::vector<double>& summed_strategy = G_CFR_STRATEGY_SUM.at(infoset_key);
            double total_summed_strategy = 0.0;
            for (double prob_sum : summed_strategy) {
                total_summed_strategy += prob_sum;
            }
            
            if (total_summed_strategy > 1e-9) {
                std::vector<double> avg_strategy(summed_strategy.size());
                for (size_t a = 0; a < summed_strategy.size(); ++a) {
                    avg_strategy[a] = summed_strategy[a] / total_summed_strategy;
                }
                profile.strategies[infoset_key] = avg_strategy;
            } else {
                // Uniform strategy if no data
                if (!actions_for_infoset.empty()) {
                    double uniform_prob = 1.0 / static_cast<double>(actions_for_infoset.size());
                    profile.strategies[infoset_key] = std::vector<double>(actions_for_infoset.size(), uniform_prob);
                }
            }
        }
    }
    
    return profile;
}

// Compute best response against a given strategy
double compute_best_response_recursive(HUState current_state, 
                                      int best_response_player,
                                      const StrategyProfile& opponent_strategy,
                                      std::map<std::string, double>& best_response_values,
                                      std::mt19937& rng) {
    
    if (current_state.game_over) {
        std::vector<double> returns = current_state.returns();
        return returns[best_response_player];
    }

    if (current_state.is_chance_node()) {
        HUState next_state = current_state;
        next_state.apply_action(Action::DEAL);
        return compute_best_response_recursive(next_state, best_response_player, opponent_strategy, best_response_values, rng);
    }

    int current_player = current_state.current_player();
    std::string infoset_key = get_infoset_key_for_state(current_state);

    std::vector<std::pair<Action, double>> legal_actions = current_state.legal_actions(
        current_state.game.pot_fraction_bet_sizes,
        current_state.game.fixed_bet_sizes_bb
    );

    if (legal_actions.empty()) {
        std::vector<double> returns = current_state.returns();
        return returns[best_response_player];
    }

    if (current_player == best_response_player) {
        // Best response player: choose action that maximizes value
        double best_value = -std::numeric_limits<double>::infinity();
        
        for (const auto& action_pair : legal_actions) {
            HUState next_state = current_state;
            next_state.apply_action(action_pair.first, action_pair.second);
            
            double action_value = compute_best_response_recursive(
                next_state, best_response_player, opponent_strategy, best_response_values, rng);
            
            if (action_value > best_value) {
                best_value = action_value;
            }
        }
        
        best_response_values[infoset_key] = best_value;
        return best_value;
        
    } else {
        // Opponent player: use their strategy
        double expected_value = 0.0;
        
        // Look up opponent's strategy for this infoset
        if (opponent_strategy.strategies.count(infoset_key) && 
            opponent_strategy.infoset_actions.count(infoset_key)) {
            
            const auto& strategy = opponent_strategy.strategies.at(infoset_key);
            
            // Ensure strategy and actions match current legal actions
            if (strategy.size() == legal_actions.size()) {
                for (size_t a = 0; a < legal_actions.size(); ++a) {
                    HUState next_state = current_state;
                    next_state.apply_action(legal_actions[a].first, legal_actions[a].second);
                    
                    double action_value = compute_best_response_recursive(
                        next_state, best_response_player, opponent_strategy, best_response_values, rng);
                    
                    expected_value += strategy[a] * action_value;
                }
            } else {
                // Fallback to uniform if strategy doesn't match
                for (const auto& action_pair : legal_actions) {
                    HUState next_state = current_state;
                    next_state.apply_action(action_pair.first, action_pair.second);
                    
                    double action_value = compute_best_response_recursive(
                        next_state, best_response_player, opponent_strategy, best_response_values, rng);
                    
                    expected_value += action_value / static_cast<double>(legal_actions.size());
                }
            }
        } else {
            // No strategy found, use uniform
            for (const auto& action_pair : legal_actions) {
                HUState next_state = current_state;
                next_state.apply_action(action_pair.first, action_pair.second);
                
                double action_value = compute_best_response_recursive(
                    next_state, best_response_player, opponent_strategy, best_response_values, rng);
                
                expected_value += action_value / static_cast<double>(legal_actions.size());
            }
        }
        
        return expected_value;
    }
}

// Compute expected utility when two strategies play against each other
double compute_strategy_vs_strategy_utility(const HUState& start_state,
                                           int target_player,
                                           const StrategyProfile& strategy_p0,
                                           const StrategyProfile& strategy_p1,
                                           int num_simulations = 1000) {
    
    double total_utility = 0.0;
    std::mt19937 rng(std::random_device{}());
    
    for (int sim = 0; sim < num_simulations; ++sim) {
        HUState current_state = start_state;
        
        while (!current_state.game_over) {
            if (current_state.is_chance_node()) {
                current_state.apply_action(Action::DEAL);
                continue;
            }
            
            int current_player = current_state.current_player();
            std::string infoset_key = get_infoset_key_for_state(current_state);
            
            std::vector<std::pair<Action, double>> legal_actions = current_state.legal_actions(
                current_state.game.pot_fraction_bet_sizes,
                current_state.game.fixed_bet_sizes_bb
            );
            
            if (legal_actions.empty()) {
                break;
            }
            
            const StrategyProfile& current_strategy = (current_player == 0) ? strategy_p0 : strategy_p1;
            
            // Sample action according to strategy
            int chosen_action_idx = 0;
            if (current_strategy.strategies.count(infoset_key)) {
                const auto& strategy = current_strategy.strategies.at(infoset_key);
                if (strategy.size() == legal_actions.size()) {
                    std::discrete_distribution<int> action_dist(strategy.begin(), strategy.end());
                    chosen_action_idx = action_dist(rng);
                } else {
                    // Fallback to uniform
                    std::uniform_int_distribution<int> uniform_dist(0, legal_actions.size() - 1);
                    chosen_action_idx = uniform_dist(rng);
                }
            } else {
                // No strategy found, use uniform
                std::uniform_int_distribution<int> uniform_dist(0, legal_actions.size() - 1);
                chosen_action_idx = uniform_dist(rng);
            }
            
            current_state.apply_action(legal_actions[chosen_action_idx].first, legal_actions[chosen_action_idx].second);
        }
        
        std::vector<double> returns = current_state.returns();
        total_utility += returns[target_player];
    }
    
    return total_utility / static_cast<double>(num_simulations);
}

// Main function to compute exploitability according to the user's specifications
double compute_exploitability(const HUState& start_state,
                             const StrategyProfile& strategy,
                             int num_simulations_per_evaluation = 1000) {
    
    std::cout << "Computing exploitability using the specified formula..." << std::endl;
    
    // Check if this is a raked game
    bool is_raked_game = start_state.game.rake_percentage > 0.0;
    std::cout << "Game type: " << (is_raked_game ? "Raked" : "Non-raked") 
              << " (rake = " << (start_state.game.rake_percentage * 100.0) << "%)" << std::endl;
    
    std::mt19937 rng(std::random_device{}());
    
    // Compute MES_EV (Maximally Exploitative Strategy Expected Value) for both players
    std::map<std::string, double> br_values_p0;
    double mes_ev_p0 = 0.0;
    
    for (int sim = 0; sim < num_simulations_per_evaluation; ++sim) {
        HUState sim_state = start_state;
        mes_ev_p0 += compute_best_response_recursive(sim_state, 0, strategy, br_values_p0, rng);
    }
    mes_ev_p0 /= static_cast<double>(num_simulations_per_evaluation);
    
    std::map<std::string, double> br_values_p1;
    double mes_ev_p1 = 0.0;
    
    for (int sim = 0; sim < num_simulations_per_evaluation; ++sim) {
        HUState sim_state = start_state;
        mes_ev_p1 += compute_best_response_recursive(sim_state, 1, strategy, br_values_p1, rng);
    }
    mes_ev_p1 /= static_cast<double>(num_simulations_per_evaluation);
    
    std::cout << "MES_EV[Player0] (Best Response EV): " << std::fixed << std::setprecision(6) << mes_ev_p0 << std::endl;
    std::cout << "MES_EV[Player1] (Best Response EV): " << std::fixed << std::setprecision(6) << mes_ev_p1 << std::endl;
    
    double exploitability;
    
    if (is_raked_game) {
        // For raked games: Exploitability = ((MES_EV[Player0] - Current_EV[Player0]) + (MES_EV[Player1] - Current_EV[Player1])) × 0.5
        
        // Compute Current_EV (expected value of the current strategy for each player)
        double current_ev_p0 = compute_strategy_vs_strategy_utility(start_state, 0, strategy, strategy, num_simulations_per_evaluation);
        double current_ev_p1 = compute_strategy_vs_strategy_utility(start_state, 1, strategy, strategy, num_simulations_per_evaluation);
        
        std::cout << "Current_EV[Player0] (Self-play EV): " << std::fixed << std::setprecision(6) << current_ev_p0 << std::endl;
        std::cout << "Current_EV[Player1] (Self-play EV): " << std::fixed << std::setprecision(6) << current_ev_p1 << std::endl;
        
        double improvement_p0 = mes_ev_p0 - current_ev_p0;
        double improvement_p1 = mes_ev_p1 - current_ev_p1;
        
        std::cout << "Potential improvement Player0: " << std::fixed << std::setprecision(6) << improvement_p0 << std::endl;
        std::cout << "Potential improvement Player1: " << std::fixed << std::setprecision(6) << improvement_p1 << std::endl;
        
        exploitability = (improvement_p0 + improvement_p1) * 0.5;
        
        std::cout << "Raked Game Exploitability = ((MES_EV[P0] - Current_EV[P0]) + (MES_EV[P1] - Current_EV[P1])) × 0.5" << std::endl;
        std::cout << "                           = ((" << mes_ev_p0 << " - " << current_ev_p0 << ") + (" 
                  << mes_ev_p1 << " - " << current_ev_p1 << ")) × 0.5" << std::endl;
        std::cout << "                           = (" << improvement_p0 << " + " << improvement_p1 << ") × 0.5" << std::endl;
        
    } else {
        // For non-raked games: Exploitability = (MES_EV[Player0] + MES_EV[Player1]) × 0.5
        exploitability = (mes_ev_p0 + mes_ev_p1) * 0.5;
        
        std::cout << "Non-raked Game Exploitability = (MES_EV[Player0] + MES_EV[Player1]) × 0.5" << std::endl;
        std::cout << "                               = (" << mes_ev_p0 << " + " << mes_ev_p1 << ") × 0.5" << std::endl;
    }
    
    std::cout << "Final Exploitability: " << std::fixed << std::setprecision(6) << exploitability << std::endl;
    
    return exploitability;
}

// Function to compute exploitability of current CFR strategy
double compute_current_strategy_exploitability(const HUState& start_state, int num_simulations) {
    StrategyProfile current_strategy = get_current_strategy_profile();
    return compute_exploitability(start_state, current_strategy, num_simulations);
}

// Function to compute exploitability from a saved strategy file
double compute_strategy_file_exploitability(const HUState& start_state, 
                                           const std::string& strategy_csv_filename,
                                           int num_simulations) {
    StrategyProfile strategy = load_strategy_from_csv(strategy_csv_filename);
    return compute_exploitability(start_state, strategy, num_simulations);
}

// ===============================
// END EXPLOITABILITY CALCULATION
// ===============================

// Outcome Sampling MCCFR - faster per iteration for large games
double mccfr_outcome_sampling(HUState current_state, std::vector<double> reach_probs, std::mt19937& rng, double sample_prob = 1.0) {
    if (current_state.game_over) {
        // Return utility for the updating player, weighted by sample probability
        std::vector<double> returns = current_state.returns();
        // In outcome sampling, we track utilities from the perspective of both players
        // Return the utility difference scaled by sample probability for proper importance sampling
        return (returns[0] - returns[1]) / sample_prob;
    }

    if (current_state.is_chance_node()) {
        HUState next_state = current_state;
        next_state.apply_action(Action::DEAL);
        return mccfr_outcome_sampling(next_state, reach_probs, rng, sample_prob);
    }

    int current_player = current_state.current_player();
    std::string infoset_key = get_infoset_key_for_state(current_state);

    std::vector<std::pair<Action, double>> legal_actions;
    int num_actions_for_infoset;

    { // Scope for G_CFR_INFOSET_ACTIONS_MUTEX
        std::lock_guard<std::mutex> lock(G_CFR_INFOSET_ACTIONS_MUTEX);
        auto it_actions = G_CFR_INFOSET_ACTIONS.find(infoset_key);
        if (it_actions == G_CFR_INFOSET_ACTIONS.end()) {
            legal_actions = current_state.legal_actions(
                current_state.game.pot_fraction_bet_sizes,
                current_state.game.fixed_bet_sizes_bb
            );
            if (legal_actions.empty()) {
                std::vector<double> returns = current_state.returns();
                return (returns[0] - returns[1]) / sample_prob;
            }
            G_CFR_INFOSET_ACTIONS[infoset_key] = legal_actions;
            G_CFR_REGRET_SUM.try_emplace(infoset_key, std::vector<double>(legal_actions.size(), 0.0));
            G_CFR_STRATEGY_SUM.try_emplace(infoset_key, std::vector<double>(legal_actions.size(), 0.0));
            num_actions_for_infoset = legal_actions.size();
        } else {
            legal_actions = it_actions->second;
            num_actions_for_infoset = legal_actions.size();
        }
    }

    if (num_actions_for_infoset == 0) { 
        std::vector<double> returns = current_state.returns();
        return (returns[0] - returns[1]) / sample_prob;
    }

    std::vector<double> strategy = get_mccfr_strategy(infoset_key, num_actions_for_infoset);

    // Sample one action according to current strategy
    std::discrete_distribution<int> action_dist(strategy.begin(), strategy.end());
    int sampled_action = action_dist(rng);
    
    HUState next_state = current_state;
    next_state.apply_action(legal_actions[sampled_action].first, legal_actions[sampled_action].second);
    
    std::vector<double> next_reach_probs = reach_probs;
    next_reach_probs[current_player] *= strategy[sampled_action];
    
    // Continue with the sampled action, updating sample probability
    double new_sample_prob = sample_prob * strategy[sampled_action];
    double utility = mccfr_outcome_sampling(next_state, next_reach_probs, rng, new_sample_prob);
    
    // The utility from current player's perspective
    double player_utility = (current_player == 0) ? utility : -utility;
    
    // Update strategy sum (weighted by reach probability of current player)
    {
        std::lock_guard<std::mutex> lock(G_CFR_STRATEGY_SUM_MUTEX);
        double weight = G_USE_LINEAR_CFR ? std::max(1, G_CFR_CURRENT_ITERATION.load()) : 1.0;
        for (int a = 0; a < num_actions_for_infoset; ++a) {
            G_CFR_STRATEGY_SUM[infoset_key][a] += weight * reach_probs[current_player] * strategy[a];
        }
    }
    
    // Update regrets using proper outcome sampling
    // In outcome sampling, regrets are updated only at nodes belonging to the updating player
    // We determine updating player by alternating or other scheme - here we'll update for both players
    {
        std::lock_guard<std::mutex> lock(G_CFR_REGRET_SUM_MUTEX);
        
        // Calculate expected utility under current strategy (baseline for regret)
        double expected_utility = player_utility; // Since we only sampled one action, this is our estimate
        
        // Weight for regret updates based on reach probabilities and sample probability
        double regret_weight = reach_probs[1 - current_player] / sample_prob;
        
        // Update regrets for each action
        for (int a = 0; a < num_actions_for_infoset; ++a) {
            double instant_regret;
            
            if (a == sampled_action) {
                // For the sampled action, regret is 0 since action_utility = expected_utility
                instant_regret = 0.0;
            } else {
                // For unsampled actions, we estimate negative regret
                // This is a key insight: in outcome sampling, unsampled actions get negative regret
                // proportional to the current utility, weighted by their probability
                instant_regret = -regret_weight * player_utility;
            }
            
            // CFR+: Clamp regrets to be non-negative
            G_CFR_REGRET_SUM[infoset_key][a] = std::max(G_CFR_REGRET_SUM[infoset_key][a] + instant_regret, 0.0);
        }
    }
    
    return utility;
}

// Accelerated MCCFR Solver with multiple optimization options
void solve_from_state_mccfr_accelerated(const HUState& start_state, 
                                       int num_iterations, 
                                       const std::string& output_csv_filename, 
                                       int num_threads_to_use,
                                       bool use_outcome_sampling,
                                       bool use_linear_cfr,
                                       int exploitability_check_interval) {
    
    // Clear global CFR data structures for a new solve
    G_CFR_REGRET_SUM.clear();
    G_CFR_STRATEGY_SUM.clear();
    G_CFR_INFOSET_ACTIONS.clear();
    G_CFR_COMPLETED_ITERATIONS_COUNT = 0;
    G_CFR_CURRENT_ITERATION = 0;
    G_USE_LINEAR_CFR = use_linear_cfr;
    
    // Reserve space for better performance
    G_CFR_REGRET_SUM.reserve(10000);
    G_CFR_STRATEGY_SUM.reserve(10000);
    G_CFR_INFOSET_ACTIONS.reserve(10000);

    std::cout << "Starting Accelerated MCCFR for " << num_iterations << " iterations using " << num_threads_to_use << " threads." << std::endl;
    std::cout << "Optimizations enabled:" << std::endl;
    std::cout << "  - CFR+ (regret clamping): YES" << std::endl;
    std::cout << "  - Linear CFR weighting: " << (use_linear_cfr ? "YES" : "NO") << std::endl;
    std::cout << "  - Sampling method: " << (use_outcome_sampling ? "Outcome Sampling" : "External Sampling") << std::endl;
    std::cout << "  - Optimized data structures: YES (unordered_map)" << std::endl;
    if (exploitability_check_interval > 0) {
        std::cout << "  - Exploitability tracking: Every " << exploitability_check_interval << " iterations" << std::endl;
    }
    std::cout << "Output will be saved to: " << output_csv_filename << std::endl;

    auto start_time = std::chrono::high_resolution_clock::now();

    std::vector<std::thread> threads;
    int iterations_per_thread_base = (num_iterations > 0 && num_threads_to_use > 0) ? (num_iterations / num_threads_to_use) : 0;
    int iterations_remainder = (num_iterations > 0 && num_threads_to_use > 0) ? (num_iterations % num_threads_to_use) : 0;

    // Exploitability tracking
    std::vector<std::tuple<int, double, double>> exploitability_history;
    std::mutex exploitability_history_mutex;
    std::atomic<int> next_exploitability_check(exploitability_check_interval);

    auto accelerated_mccfr_task = [&](int iterations_for_this_thread, int thread_id) {
        std::mt19937 thread_rng(std::random_device{}() + thread_id);
        
        for (int i = 0; i < iterations_for_this_thread; ++i) {
            // Increment iteration counter for Linear CFR weighting
            int current_iter = G_CFR_CURRENT_ITERATION.fetch_add(1, std::memory_order_relaxed);
            
            if (use_outcome_sampling) {
                // Use outcome sampling MCCFR
                std::vector<double> initial_reach_probs(2, 1.0);
                mccfr_outcome_sampling(start_state, initial_reach_probs, thread_rng);
            } else {
                // Use external sampling MCCFR
                int updating_player = current_iter % 2;
                mccfr_recursive(start_state, updating_player, thread_rng);
            }
            
            int completed_count = G_CFR_COMPLETED_ITERATIONS_COUNT.fetch_add(1, std::memory_order_relaxed) + 1;
            
            // Check if this thread should perform an exploitability check
            if (exploitability_check_interval > 0) {
                int expected_check = next_exploitability_check.load(std::memory_order_acquire);
                if (completed_count >= expected_check && completed_count > 0) {
                    if (next_exploitability_check.compare_exchange_strong(expected_check, expected_check + exploitability_check_interval)) {
                        std::cout << "\n--- Accelerated MCCFR Exploitability Check at Iteration " << completed_count << " ---" << std::endl;
                        
                        try {
                            double exploitability = compute_current_strategy_exploitability(start_state, 500);
                            double starting_pot = start_state.cumulative_pot;
                            double exploitability_percentage = (starting_pot > 0) ? (exploitability / starting_pot) * 100.0 : 0.0;
                            
                            {
                                std::lock_guard<std::mutex> lock(exploitability_history_mutex);
                                exploitability_history.push_back(std::make_tuple(completed_count, exploitability, exploitability_percentage));
                            }
                            
                            std::cout << "Iteration " << completed_count << " - Exploitability: " << std::fixed 
                                      << std::setprecision(6) << exploitability << " BB/hand (" 
                                      << std::fixed << std::setprecision(2) << exploitability_percentage << "% of pot)" << std::endl;
                            
                        } catch (const std::exception& e) {
                            std::cerr << "Error computing exploitability: " << e.what() << std::endl;
                        }
                    }
                }
            }
        }
    };

    for (int t = 0; t < num_threads_to_use; ++t) {
        int iterations_for_this_thread = iterations_per_thread_base + (t < iterations_remainder ? 1 : 0);
        if (iterations_for_this_thread > 0) {
            threads.emplace_back(accelerated_mccfr_task, iterations_for_this_thread, t);
        }
    }

    // Progress monitoring
    int report_interval_ms = 1000; // Report every second
    int iterations_reported_at_last_print = -1;
    bool first_report_triggered = false;

    if (num_iterations > 0) {
        std::cout << "Accelerated MCCFR progress (target " << num_iterations << "):" << std::endl;
    }

    while(G_CFR_COMPLETED_ITERATIONS_COUNT.load(std::memory_order_acquire) < num_iterations) {
        std::this_thread::sleep_for(std::chrono::milliseconds(report_interval_ms));
        int current_completed = G_CFR_COMPLETED_ITERATIONS_COUNT.load(std::memory_order_acquire);
        
        if (!first_report_triggered || current_completed > iterations_reported_at_last_print) {
            auto current_time = std::chrono::high_resolution_clock::now();
            auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(current_time - start_time);
            double iterations_per_second = (elapsed.count() > 0) ? (current_completed * 1000.0 / elapsed.count()) : 0.0;
            
            std::cout << "  Completed " << current_completed << " / " << num_iterations 
                      << " iterations (" << std::fixed << std::setprecision(1) << iterations_per_second 
                      << " iter/sec)" << std::endl;
            iterations_reported_at_last_print = current_completed;
            first_report_triggered = true;
        }
        if (current_completed >= num_iterations) break; 
    }

    // Join threads
    for (auto& th : threads) {
        if (th.joinable()) {
            th.join();
        }
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto total_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    double final_iterations_per_second = (total_elapsed.count() > 0) ? (num_iterations * 1000.0 / total_elapsed.count()) : 0.0;
    
    std::cout << "Accelerated MCCFR complete in " << total_elapsed.count() << "ms" << std::endl;
    std::cout << "Average performance: " << std::fixed << std::setprecision(1) 
              << final_iterations_per_second << " iterations/second" << std::endl;

    // Save strategy with performance info
    std::cout << "Saving accelerated strategy..." << std::endl;
    std::ofstream outfile(output_csv_filename);
    outfile << "infoset_key,action_string,probability\n";

    for (const auto& pair : G_CFR_INFOSET_ACTIONS) {
        const std::string& infoset_key = pair.first;
        const auto& actions_for_infoset = pair.second;

        if (G_CFR_STRATEGY_SUM.count(infoset_key)) {
            const std::vector<double>& summed_strategy = G_CFR_STRATEGY_SUM.at(infoset_key);
            double total_summed_strategy = 0.0;
            for (double prob_sum : summed_strategy) {
                total_summed_strategy += prob_sum;
            }

            if (total_summed_strategy > 1e-9) {
                for (size_t a = 0; a < summed_strategy.size(); ++a) {
                    double avg_prob = summed_strategy[a] / total_summed_strategy;
                    outfile << "\"" << infoset_key << "\",\"" 
                            << action_to_string(actions_for_infoset[a].first, actions_for_infoset[a].second) 
                            << "\"," << std::fixed << std::setprecision(5) << avg_prob << "\n";
                }
            } else {
                if (!actions_for_infoset.empty()) {
                    double uniform_prob = 1.0 / static_cast<double>(actions_for_infoset.size());
                    for (size_t a = 0; a < actions_for_infoset.size(); ++a) {
                        outfile << "\"" << infoset_key << "\",\""
                                << action_to_string(actions_for_infoset[a].first, actions_for_infoset[a].second)
                                << "\"," << std::fixed << std::setprecision(5) << uniform_prob << "\n";
                    }
                }
            }
        }
    }
    outfile.close();
    
    std::cout << "Accelerated strategy saved to " << output_csv_filename << std::endl;
    std::cout << "Found " << G_CFR_INFOSET_ACTIONS.size() << " unique information sets" << std::endl;
}
