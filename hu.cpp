#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <random>
#include <chrono>
#include <sstream>
#include <set>
#include <stdexcept>
#include <numeric>
#include <iterator>
#include <unordered_map>
#include <map>
#include <cassert>
#include <utility> // For std::pair
#include <iomanip> // for setprecision
#include <functional> // for std::function and std::greater
#include <tuple> // For std::tuple in round_action_history

// Constants definitions
// constexpr int NUM_PLAYERS = 2;
// constexpr double INITIAL_STACK = 100.0; // Example: 100 Big Blinds
// constexpr double SMALL_BLIND_AMOUNT = 0.5;
// constexpr double BIG_BLIND_AMOUNT = 1.0;
// static const float RAKE_PERCENTAGE = 0.07;  // 7% rake (currently not applied in returns)

//----------------------------------------------------------------------------
// Card, Deck, and Helper Functions
//----------------------------------------------------------------------------

struct Card {
    std::string rank;
    std::string suit;

    Card() = default;
    
    Card(const std::string& r, char s) : rank(r), suit(1, s) {}

    std::string toString() const {
        return rank + suit;
    }

    // For sorting and set operations if needed
    bool operator<(const Card& other) const {
        if (rank != other.rank) return rank < other.rank; // Simplistic comparison, can be improved
        return suit < other.suit;
    }
};

const std::vector<std::string> RANKS = {"2", "3", "4", "5", "6", "7", "8", "9", "10", "J", "Q", "K", "A"};
const std::vector<char> SUITS_CHARS = {'h', 'd', 'c', 's'}; // Renamed to avoid conflict if SUITS is used elsewhere

std::vector<Card> make_deck() {
    std::vector<Card> deck;
    for (const auto& rank : RANKS) {
        for (const auto& suit_char : SUITS_CHARS) {
            Card card;
            card.rank = rank;
            card.suit = std::string(1, suit_char);
            deck.push_back(card);
        }
    }
    return deck;
}

// ----------------------------------------------------------------------------
// PokerEvaluator (largely similar to spingo.cpp)
// ----------------------------------------------------------------------------
class PokerEvaluator {
public:
    const std::unordered_map<std::string, int> RANK_VALUES = {
        {"2", 2}, {"3", 3}, {"4", 4}, {"5", 5},
        {"6", 6}, {"7", 7}, {"8", 8}, {"9", 9},
        {"10", 10}, {"J", 11}, {"Q", 12}, {"K", 13}, {"A", 14}
    };

    enum HandRank {
        HIGH_CARD = 0, PAIR = 1, TWO_PAIR = 2, THREE_KIND = 3,
        STRAIGHT = 4, FLUSH = 5, FULL_HOUSE = 6, FOUR_KIND = 7,
        STRAIGHT_FLUSH = 8
    };

    int evaluateHand(const std::vector<Card>& holeCards, const std::vector<Card>& communityCards) {
        std::vector<Card> allCards = holeCards;
        allCards.insert(allCards.end(), communityCards.begin(), communityCards.end());
        
        if (allCards.size() < 5) return 0; // Not enough cards to form a hand

        std::vector<std::vector<Card>> combinations;
        std::vector<Card> currentCombo;
        generateCombinations(allCards, 5, 0, currentCombo, combinations);
        
        int bestScore = 0;
        for (const auto &combo : combinations) {
            int score = evaluateFiveCardHand(combo);
            if (score > bestScore) {
                bestScore = score;
            }
        }
        return bestScore;
    }

    int evaluateFiveCardHand(const std::vector<Card>& hand) {
        std::vector<int> cardValues;
        for (const auto &c : hand)
            cardValues.push_back(RANK_VALUES.at(c.rank));
        std::sort(cardValues.rbegin(), cardValues.rend());

        std::map<int, int, std::greater<int>> freq;
        for (int v : cardValues)
            freq[v]++;

        bool isFlush = checkFlush(hand);
        int straightHigh = 0;
        bool isStraight = checkStraight(hand, straightHigh);

        if (isFlush && isStraight) return makeScore(HandRank::STRAIGHT_FLUSH, {straightHigh});
        
        for (const auto &entry : freq) {
            if (entry.second == 4) {
                int kicker = 0;
                for (int v : cardValues) {
                    if (v != entry.first) { kicker = v; break; }
                }
                return makeScore(HandRank::FOUR_KIND, {entry.first, kicker});
            }
        }
        
        int tripleRank = 0, pairRank = 0;
        for (const auto &entry : freq) {
            if (entry.second >= 3 && tripleRank == 0)
                tripleRank = entry.first;
        }
        if (tripleRank) {
            for (const auto &entry : freq) {
                if (entry.first != tripleRank && entry.second >= 2) {
                    pairRank = entry.first;
                    break;
                }
            }
            if (pairRank)
                return makeScore(HandRank::FULL_HOUSE, {tripleRank, pairRank});
        }
        
        if (isFlush) return makeScore(HandRank::FLUSH, cardValues);
        if (isStraight) return makeScore(HandRank::STRAIGHT, {straightHigh});
        
        for (const auto &entry : freq) {
            if (entry.second == 3) {
                std::vector<int> kickers;
                for (int v : cardValues) {
                    if (v != entry.first) kickers.push_back(v);
                }
                return makeScore(HandRank::THREE_KIND, {entry.first, kickers[0], kickers[1]});
            }
        }
        
        std::vector<int> pairs;
        int kicker = 0;
        for (const auto &entry : freq) {
            if (entry.second >= 2)
                pairs.push_back(entry.first);
        }
        if (pairs.size() >= 2) {
            // Already sorted by rank due to map properties, take top 2
            for (int v : cardValues) {
                if (v != pairs[0] && v != pairs[1]) { kicker = v; break; }
            }
            return makeScore(HandRank::TWO_PAIR, {pairs[0], pairs[1], kicker});
        }
        
        for (const auto &entry : freq) {
            if (entry.second == 2) {
                std::vector<int> kickers;
                for (int v : cardValues) {
                    if (v != entry.first) kickers.push_back(v);
                }
                return makeScore(HandRank::PAIR, {entry.first, kickers[0], kickers[1], kickers[2]});
            }
        }
        
        return makeScore(HandRank::HIGH_CARD, cardValues);
    }

    int makeScore(int handRank, const std::vector<int>& vals) {
        long long score = static_cast<long long>(handRank) * 10000000000LL; // Ensure enough space
        long long multiplier = 100000000;
        for (int v : vals) {
            if (multiplier < 1) break; // Max 5 kickers
            score += static_cast<long long>(v) * multiplier;
            multiplier /= 100;
        }
        // Since original scores were int, let's cap or check for overflow if necessary,
        // but for ranking, long long is safer. If int is required, scale down.
        // For simplicity, returning as int, assuming values fit.
        if (score > std::numeric_limits<int>::max()) return std::numeric_limits<int>::max();
        return static_cast<int>(score);
    }

    void generateCombinations(const std::vector<Card>& cards,
                              size_t combinationSize,
                              size_t start,
                              std::vector<Card>& currentCombo,
                              std::vector<std::vector<Card>>& allCombinations) {
        if (currentCombo.size() == combinationSize) {
            allCombinations.push_back(currentCombo);
            return;
        }
        size_t needed = combinationSize - currentCombo.size();
        for (size_t i = start; i < cards.size() && (cards.size() - i >= needed) ; ++i) {
            currentCombo.push_back(cards[i]);
            generateCombinations(cards, combinationSize, i + 1, currentCombo, allCombinations);
            currentCombo.pop_back();
        }
    }

    bool checkFlush(const std::vector<Card>& hand) {
        if (hand.empty()) return false;
        std::string suit = hand.front().suit;
        return std::all_of(hand.begin(), hand.end(), [&](const Card &c) { return c.suit == suit; });
    }

    bool checkStraight(const std::vector<Card>& hand, int &highStraightValue) {
        if (hand.size() != 5) return false;
        std::vector<int> values;
        for (const auto &card : hand)
            values.push_back(RANK_VALUES.at(card.rank));
        std::sort(values.begin(), values.end());
        
        // Remove duplicates for straight check
        values.erase(std::unique(values.begin(), values.end()), values.end());
        if (values.size() < 5) return false; // Not enough unique ranks for a 5-card straight

        bool normalStraight = true;
        for (size_t i = 0; i < 4; ++i) {
            if (values[i+1] - values[i] != 1) {
                normalStraight = false;
                break;
            }
        }
        if (normalStraight) {
            highStraightValue = values[4];
            return true;
        }
        
        // Check for A-5 straight (wheel)
        bool isWheel = true;
        std::vector<int> wheel_ranks = {14, 2, 3, 4, 5}; // A, 2, 3, 4, 5
        std::vector<int> sorted_hand_ranks = values; // Already sorted
        if (sorted_hand_ranks.size() == 5 &&
            sorted_hand_ranks[0] == RANK_VALUES.at("2") &&
            sorted_hand_ranks[1] == RANK_VALUES.at("3") &&
            sorted_hand_ranks[2] == RANK_VALUES.at("4") &&
            sorted_hand_ranks[3] == RANK_VALUES.at("5") &&
            sorted_hand_ranks[4] == RANK_VALUES.at("A")) {
            highStraightValue = 5; // Ace plays low
            return true;
        }
        return false;
    }
};

// ----------------------------------------------------------------------------
// Action enumeration and action amounts mapping
// ----------------------------------------------------------------------------

enum class Action {
    UNKNOWN = -1,
    DEAL = 0,       // Chance action
    FOLD = 1,
    CHECK = 2,      // CHECK is also CALL when current_bet_to_match is 0 and player has 0 in pot.
    CALL = 3,
    BET_RAISE = 4,  // This will now be associated with a specific amount from legal_actions
    ALL_IN = 5,
    POST_SB = 6,    // Internal action for posting small blind
    POST_BB = 7     // Internal action for posting big blind
};




std::string action_to_string(Action action, double amount = 0.0) {
    switch (action) {
        case Action::DEAL: return "DEAL";
        case Action::FOLD: return "FOLD";
        case Action::CHECK: return "CHECK";
        case Action::CALL: return "CALL";
        case Action::BET_RAISE: {
            std::ostringstream oss;
            oss << "BET_RAISE(" << std::fixed << std::setprecision(2) << amount << ")";
            return oss.str();
        }
        case Action::ALL_IN: return "ALL_IN";
        case Action::POST_SB: return "POST_SB";
        case Action::POST_BB: return "POST_BB";
        default: return "UNKNOWN";
    }
}

// ----------------------------------------------------------------------------
// HUGame class (Heads-Up Game)
// ----------------------------------------------------------------------------
class HUState;

class HUGame {
public:
    const int num_players;
    const double initial_stack;
    const double small_blind_amount;
    const double big_blind_amount;
    const float rake_percentage;
    const std::vector<double> pot_fraction_bet_sizes;
    const std::vector<double> fixed_bet_sizes_bb;
    const double all_in_threshold;        // Default 150% (1.5)
    const double force_all_in_threshold;  // Default 20% (0.2)
    const double merging_threshold;       // Default 10% (0.1)

    HUGame(
        int num_players_val = 2,
        double initial_stack_val = 100.0,
        double small_blind_amount_val = 0.5,
        double big_blind_amount_val = 1.0,
        float rake_percentage_val = 0.07,
        const std::vector<double>& pot_fractions = {0.33, 0.5, 0.75, 1.0},
        const std::vector<double>& fixed_bbs = {2.5, 3.0, 4.0},
        double all_in_threshold_val = 1.5,
        double force_all_in_threshold_val = 0.2,
        double merging_threshold_val = 0.1
    );

    HUState new_initial_state() const;
};

// ----------------------------------------------------------------------------
// HUState class (Heads-Up State)
// ----------------------------------------------------------------------------

class HUState {
public:
    const HUGame& game;  // Reference to the game configuration

    // Add action_fixed_amounts to member variables
    std::unordered_map<Action, double> action_fixed_amounts;

    // New attributes
    std::vector<std::pair<Action, double>> bets;  // Replace current_round_bets
    double pot;  // Current pot for the round
    double cumulative_pot;  // Total pot across all rounds
    std::vector<int> active_players;  // Replace active_players_indices set
    double current_bet;  // Replace current_bet_to_match

    // Keep old attributes during transition
    std::vector<Action> current_round_bets;  // Will be replaced by bets
    std::set<int> active_players_indices;    // Will be replaced by active_players
    double current_bet_to_match;             // Will be replaced by current_bet

    // Existing attributes to keep
    std::string round;
    std::map<std::string, std::vector<std::tuple<int, Action, double>>> round_action_history;
    std::vector<Card> cards;
    std::vector<double> pot_contribution_this_round;
    std::vector<double> players_stack;
    bool game_over;
    int next_player_idx;
    std::vector<Card> community_cards;
    std::vector<double> cumulative_pot_contribution;
    std::vector<Card> deck;
    double last_bet_or_raise_increment_this_round;
    int last_aggressor_idx_this_round;
    std::mt19937 rng;

    HUState(const HUGame& game_ref) : game(game_ref) {
        // Initialize action_fixed_amounts
        action_fixed_amounts = {
            {Action::POST_SB, game.small_blind_amount},
            {Action::POST_BB, game.big_blind_amount}
        };

        // Initialize new attributes
        pot = game.small_blind_amount + game.big_blind_amount;
        cumulative_pot = game.small_blind_amount + game.big_blind_amount;
        active_players = {0, 1};
        current_bet = game.big_blind_amount;

        // Update ACTION_FIXED_AMOUNTS to use game parameters
        std::unordered_map<Action, double> action_fixed_amounts = {
            {Action::POST_SB, game.small_blind_amount},
            {Action::POST_BB, game.big_blind_amount}
        };
        // Initialize old attributes (during transition)
        current_round_bets = std::vector<Action>(game.num_players, Action::UNKNOWN);
        active_players_indices = std::set<int>{0, 1};
        current_bet_to_match = current_bet;

        // Rest of initialization using game parameters
        pot_contribution_this_round = std::vector<double>(game.num_players, 0.0);
        players_stack = std::vector<double>(game.num_players, game.initial_stack);
        cumulative_pot_contribution = std::vector<double>(game.num_players, 0.0);
        game_over = false;

        // Post SB
        pot_contribution_this_round[0] = game.small_blind_amount;
        players_stack[0] -= game.small_blind_amount;
        current_round_bets[0] = Action::POST_SB;

        // Post BB
        pot_contribution_this_round[1] = game.big_blind_amount;
        players_stack[1] -= game.big_blind_amount;
        current_round_bets[1] = Action::POST_BB;

        next_player_idx = 0;
        last_bet_or_raise_increment_this_round = game.big_blind_amount;
        last_aggressor_idx_this_round = 1;

        deck = make_deck();
        rng.seed(std::random_device{}());
        std::shuffle(deck.begin(), deck.end(), rng);

        round = "preflop";
    }

    double get_fixed_amount(Action a) {
        if(action_fixed_amounts.count(a))
            return action_fixed_amounts[a];
        return 0.0;
    }
    
    // Helper method to keep old and new attributes in sync
    void sync_attributes() {
        // Sync active players
        active_players_indices = std::set<int>(active_players.begin(), active_players.end());
        
        // Sync bets (no longer needed as current_round_bets is canonical)
        // for (int i = 0; i < game.num_players; i++) {
        //     current_round_bets[i] = bets[i].first;
        // }
        
        // Sync current bet
        current_bet = current_bet_to_match;
    }

    int current_player() const {
        if (game_over) return -1;
        return next_player_idx;
    }

    bool is_chance_node() const {
        if (round == "preflop" && cards.size() < game.num_players * 2) return true;
        if (round == "flop" && community_cards.size() < 3) return true;
        if (round == "turn" && community_cards.size() < 4) return true;
        if (round == "river" && community_cards.size() < 5) return true;
        return false;
    }

    void deal_cards() {
        if (round == "preflop" && cards.empty()) {
            // Deal hole cards
            for (int p = 0; p < game.num_players; ++p) {
                cards.push_back(deck.back()); deck.pop_back();
                cards.push_back(deck.back()); deck.pop_back();
            }
            next_player_idx = 0; // SB (player 0) acts first preflop
        } else if (round == "flop" && community_cards.empty()) {
            // Burn one
            if (!deck.empty()) deck.pop_back();
            // Deal 3 for flop
            for (int i = 0; i < 3; ++i) {
                if (!deck.empty()) { community_cards.push_back(deck.back()); deck.pop_back(); }
            }
            next_player_idx = (active_players_indices.count(0)) ? 0 : 1; // SB acts first if in
        } else if (round == "turn" && community_cards.size() == 3) {
            // Burn one
            if (!deck.empty()) deck.pop_back();
            // Deal 1 for turn
            if (!deck.empty()) { community_cards.push_back(deck.back()); deck.pop_back(); }
            next_player_idx = (active_players_indices.count(0)) ? 0 : 1; // SB acts first if in
        } else if (round == "river" && community_cards.size() == 4) {
            // Burn one
            if (!deck.empty()) deck.pop_back();
            // Deal 1 for river
            if (!deck.empty()) { community_cards.push_back(deck.back()); deck.pop_back(); }
            next_player_idx = (active_players_indices.count(0)) ? 0 : 1; // SB acts first if in
        }
        // If advancing to a round where no cards are dealt (e.g. showdown), next_player_idx is handled by advance_round
    }

    // Get min total amount for a bet or raise for the current player
    double get_min_bet_raise_total_amount() const {
        int player_idx = current_player();
        double current_player_contribution = pot_contribution_this_round[player_idx];
        
        // If no bet to match (i.e., first bet), minimum is BB
        if (current_bet_to_match == 0.0) {
            return current_player_contribution + game.big_blind_amount;
        }
        
        // For raises, minimum is last bet plus the size of the last bet/raise
        double min_raise_to = current_bet_to_match + last_bet_or_raise_increment_this_round;
        
        // If player already has some money in the pot this round, subtract it
        return min_raise_to;
    }
    
    // Get max bet/raise amount (which is all-in)
    double get_max_bet_raise_total_amount() const {
        int player_idx = current_player();
        return pot_contribution_this_round[player_idx] + players_stack[player_idx];
    }


    std::vector<std::pair<Action, double>> legal_actions(
        const std::vector<double>& pot_fraction_bet_sizes,
        const std::vector<double>& fixed_bet_sizes_bb
    ) {
        std::vector<std::pair<Action, double>> legal;
        if (game_over) return legal;
        if (is_chance_node()) {
            legal.push_back({Action::DEAL, 0.0});
            return legal;
        }

        int player_idx = current_player();

        // Player validity check
        if (player_idx == -1 || !active_players_indices.count(player_idx)) {
            return legal; // No current player or player not active
        }

        double player_stack = players_stack[player_idx];
        double current_player_pot_contrib = pot_contribution_this_round[player_idx];
        double amount_to_call = current_bet_to_match - current_player_pot_contrib;

        // If player has no stack and is facing a bet, they have no actions.
        if (player_stack <= 1e-5 && amount_to_call > 1e-5) {
            return legal;
        }

        bool can_check = amount_to_call <= 1e-5;

        // CHECK or FOLD (if checking)
        if (can_check) {
            legal.push_back({Action::CHECK, 0.0});
            // If current_bet_to_match is 0, FOLD is not a valid action.
            // A player can only fold if they are facing a bet.
            if (current_bet_to_match > 1e-5 && player_stack > 1e-5) {
                legal.push_back({Action::FOLD, 0.0});
            }
        } else { // Must call an amount > 0, or fold
            // FOLD is always an option if facing a bet and has stack
            // (If stack is 0 and facing bet, already returned empty)
            if (player_stack > 1e-5) {
                 legal.push_back({Action::FOLD, 0.0});
            }

            // CALL (only if not an all-in call and player has stack)
            if (player_stack > 1e-5 && amount_to_call < player_stack - 1e-5) {
                legal.push_back({Action::CALL, amount_to_call});
            }
            // If amount_to_call >= player_stack, it's an all-in call, handled by ALL_IN action.
        }

        // BET_RAISE and ALL_IN (only if player has stack)
        if (player_stack > 1e-5) {
            double all_in_total_contribution = current_player_pot_contrib + player_stack;

            // Minimum total contribution for a valid (non-all-in) bet/raise:
            double min_valid_bet_raise_total_target;
            if (current_bet_to_match == 0) { // Opening bet
                min_valid_bet_raise_total_target = current_player_pot_contrib + game.big_blind_amount;
            } else { // Raising
                min_valid_bet_raise_total_target = current_bet_to_match + std::max(game.big_blind_amount, last_bet_or_raise_increment_this_round);
            }

            std::set<double> added_bet_raise_total_amounts; // Track total contribution amounts for BET_RAISE

            // Add ALL_IN action. Amount is player's remaining stack.
            // This covers all-in calls, all-in opening bets, and all-in raises.
            legal.push_back({Action::ALL_IN, player_stack});

            // BET_RAISE options (non-all-in)
            // These are valid if they meet min raise requirements and are strictly less than all-in.
            if (all_in_total_contribution > min_valid_bet_raise_total_target + 1e-5) { // Check if a non-all-in raise is possible at all

                // Configured Pot Fraction Bets
                double current_total_pot_for_sizing = get_total_pot(); // Pot before current player acts
                for (double fraction : pot_fraction_bet_sizes) {
                    double bet_value_to_add; // Amount to add on top of current_player_pot_contrib or current_bet_to_match

                    if (current_bet_to_match == 0) { // Opening bet: fraction of current pot
                        bet_value_to_add = current_total_pot_for_sizing * fraction;
                        bet_value_to_add = std::max(bet_value_to_add, game.big_blind_amount);
                    } else { // Raising: fraction of pot *after* hypothetical call
                        double pot_if_called = current_total_pot_for_sizing + std::max(0.0, amount_to_call);
                        bet_value_to_add = pot_if_called * fraction;
                        bet_value_to_add = std::max(bet_value_to_add, std::max(game.big_blind_amount, last_bet_or_raise_increment_this_round));
                    }

                    double total_contribution_target;
                    if (current_bet_to_match == 0) {
                        total_contribution_target = current_player_pot_contrib + bet_value_to_add;
                    } else {
                        total_contribution_target = current_bet_to_match + bet_value_to_add;
                    }
                    
                    total_contribution_target = std::round(total_contribution_target / game.small_blind_amount) * game.small_blind_amount;

                    if (total_contribution_target >= min_valid_bet_raise_total_target - 1e-5 && // meet min raise/bet
                        total_contribution_target < all_in_total_contribution - 1e-5 &&      // strictly less than all-in
                        (total_contribution_target - current_player_pot_contrib <= player_stack + 1e-5) && // player has enough
                        added_bet_raise_total_amounts.find(total_contribution_target) == added_bet_raise_total_amounts.end()) {
                        
                        double amount_player_actually_adds = total_contribution_target - current_player_pot_contrib;
                        if (amount_player_actually_adds > 1e-5) { // Must add a positive amount
                             legal.push_back({Action::BET_RAISE, total_contribution_target});
                             added_bet_raise_total_amounts.insert(total_contribution_target);
                        }
                    }
                }

                // Configured Fixed BB Bets (as total bet size for the round)
                for (double bb_multiple : fixed_bet_sizes_bb) {
                    double total_contribution_target = bb_multiple * game.big_blind_amount;
                    total_contribution_target = std::round(total_contribution_target / game.small_blind_amount) * game.small_blind_amount;

                     if (total_contribution_target >= min_valid_bet_raise_total_target - 1e-5 && // meet min raise/bet
                        total_contribution_target < all_in_total_contribution - 1e-5 &&       // strictly less than all-in
                        (total_contribution_target - current_player_pot_contrib <= player_stack + 1e-5) && // player has enough
                        added_bet_raise_total_amounts.find(total_contribution_target) == added_bet_raise_total_amounts.end()) {
                        
                        double amount_player_actually_adds = total_contribution_target - current_player_pot_contrib;
                         if (amount_player_actually_adds > 1e-5) { // Must add a positive amount
                            legal.push_back({Action::BET_RAISE, total_contribution_target});
                            added_bet_raise_total_amounts.insert(total_contribution_target);
                        }
                    }
                }
            }
        }
        
        // Ensure FOLD is not present if only CHECK is possible and stack is 0 (all-in and checked to)
        if (can_check && player_stack <= 1e-5) {
            legal.erase(std::remove_if(legal.begin(), legal.end(), [](const std::pair<Action, double>& p) {
                return p.first == Action::FOLD;
            }), legal.end());
        }

        // Remove duplicate actions (e.g. if CHECK and FOLD are the only options)
        // Sort and unique can be helpful for final cleanup if complex scenarios lead to duplicates,
        // but the current logic aims to add distinct categories.
        // A simple sort and unique based on {Action, Amount} might be good.
        std::sort(legal.begin(), legal.end(), [](const auto& a, const auto& b){
            if (a.first != b.first) return a.first < b.first;
            return a.second < b.second;
        });
        legal.erase(std::unique(legal.begin(), legal.end(), [](const auto& a, const auto& b){
            return a.first == b.first && std::abs(a.second - b.second) < 1e-5;
        }), legal.end());

        return legal;
    }

    bool betting_round_complete() {
        if (active_players_indices.size() <= 1) {
            return true; // Game over if only one or zero active players
        }

        // Special case for PREFLOP when BB has the option to check/raise
        // This must be checked BEFORE the general contribution check
        if (round == "preflop" && current_bet_to_match == game.big_blind_amount) {
            // If SB has called/matched BB amount, BB still needs to act (unless BB has already acted)
            if (pot_contribution_this_round[0] == game.big_blind_amount) { // SB called or raised to BB
                // If it's BB's turn, they haven't acted yet, so round is NOT complete
                if (next_player_idx == 1) {
                    return false; // BB needs to act on their option
                }
                // If BB has explicitly checked their option, round is complete
                if (current_round_bets[1] == Action::CHECK) {
                    return true; // Round is complete after BB checks option
                }
                // If BB did something other than check (like raise), continue with general logic
            }
        }

        // Check if all active, non-all-in players have contributed equally to the current bet.
        // If someone still needs to match the current bet, the round is NOT complete.
        for (int p_idx : active_players_indices) {
            if (players_stack[p_idx] > 0 && pot_contribution_this_round[p_idx] < current_bet_to_match) {
                return false; // Someone still needs to act on the current bet.
            }
        }

        // General case for all rounds: If current_bet_to_match is 0 (all players checked or everyone folded to zero bet)
        // AND all active players have *explicitly checked* in this round.
        // This also covers post-flop two checks scenario.
        if (current_bet_to_match < 1e-5) { // Effectively 0
            bool all_active_checked_this_round = true;
            for (int p_idx : active_players_indices) {
                // If a player is active and their current action is NOT CHECK, or it's UNKNOWN (meaning they haven't acted yet in this 0-bet scenario)
                // then not all have checked.
                if (current_round_bets[p_idx] != Action::CHECK) { // Check for explicit CHECK action
                    all_active_checked_this_round = false;
                    break;
                }
            }
            if (all_active_checked_this_round && active_players_indices.size() > 0) { // All active players have explicitly checked
                return true;
            }
        }
        
        // Final fallback: if there was an aggressor, and it's their turn again, and all others matched.
        // This implicitly assumes the non-aggressor has acted and called/folded.
        // The `pot_contribution_this_round[p_idx] < current_bet_to_match` check handles if someone still needs to act.
        if (last_aggressor_idx_this_round != -1) { // There was a bet or raise this round
            // If it's the aggressor's turn again, it means the other player acted (called or folded).
            if (next_player_idx == last_aggressor_idx_this_round) {
                return true;
            }
        }

        // Otherwise, betting continues
        sync_attributes(); // Make sure attributes are in sync before returning false. (This sync is also done at start of apply_action if needed)
        return false;
    }


    void advance_to_next_player() {
        if (active_players_indices.empty()) {
            game_over = true;
            return;
        }
        // In HU, it's just the other active player.
        for (int p_idx : active_players_indices) {
            if (p_idx != next_player_idx) {
                next_player_idx = p_idx;
                return;
            }
        }
        // If only one active player, next_player_idx remains, game should end.
        if (active_players_indices.size() == 1) {
             next_player_idx = *active_players_indices.begin(); // Should be caught by game_over logic
        }
    }

    void advance_round() {
        for (int p = 0; p < game.num_players; p++) {
            cumulative_pot_contribution[p] += pot_contribution_this_round[p];
            pot_contribution_this_round[p] = 0.0;
            current_round_bets[p] = Action::UNKNOWN; // Reset actions for the new round
        }
        
        current_bet_to_match = 0.0;
        last_bet_or_raise_increment_this_round = 0.0; // BB for preflop, 0 for postflop
        last_aggressor_idx_this_round = -1;


        if (round == "preflop") {
            round = "flop";
            last_bet_or_raise_increment_this_round = 0; // No bet carried over
        } else if (round == "flop") {
            round = "turn";
        } else if (round == "turn") {
            round = "river";
        } else if (round == "river") {
            round = "showdown";
        }

        if (round != "showdown") {
            // Determine next player for post-flop rounds (SB/Button acts first if active)
            if (active_players_indices.count(0)) { // Player 0 (SB) is active
                next_player_idx = 0;
            } else if (active_players_indices.count(1)) { // Player 1 (BB) is active
                next_player_idx = 1;
            } else {
                game_over = true; // No active players
            }
        }
    }

    bool should_end_game() {
        if (active_players_indices.size() <= 1 || round == "showdown") {
            return true;
        }
        // If all remaining active players are all-in.
        bool all_active_are_all_in = true;
        if (active_players_indices.empty()) return true; // Should be caught by size <=1

        for (int p_idx : active_players_indices) {
            if (players_stack[p_idx] > 0) {
                all_active_are_all_in = false;
                break;
            }
        }
        return all_active_are_all_in;
    }

    void apply_action(Action action, double chosen_amount_for_bet_raise = 0.0) {
        double amount_to_put_in_pot = 0.0;
        
        if (is_chance_node()) {
            if (action == Action::DEAL) {
                deal_cards();
            } else {
                throw std::runtime_error("Invalid action on a chance node. Expected DEAL.");
            }
            // After dealing, check if betting round should start or if more cards needed
            if (!is_chance_node() && betting_round_complete()) { // e.g. preflop deal done, SB/BB posted
                 // This case should be handled by game flow: deal -> player actions
            }
            return;
        }

        int player_idx = current_player();
        if (player_idx == -1) throw std::runtime_error("Apply action called when game is over or no current player.");
        
        current_round_bets[player_idx] = action;
        double history_amount = 0.0; // For round_action_history

        if (action == Action::FOLD) {
            active_players_indices.erase(player_idx);
            if (active_players_indices.size() <= 1) game_over = true;
        } else if (action == Action::CHECK) {
            // No change to pot or stack, current_bet_to_match remains.
            // If this check closes the action (e.g. BB checks option, or checks around post-flop)
            // last_aggressor_idx_this_round remains -1 or the previous aggressor.
        } else if (action == Action::CALL) {
            amount_to_put_in_pot = current_bet_to_match - pot_contribution_this_round[player_idx];
            if (amount_to_put_in_pot >= players_stack[player_idx]) { // Call is all-in
                amount_to_put_in_pot = players_stack[player_idx];
                current_round_bets[player_idx] = Action::ALL_IN; // Treat as all-in
                action = Action::ALL_IN; // Update action for history
            }
            history_amount = amount_to_put_in_pot; // This is the amount added
            pot_contribution_this_round[player_idx] += amount_to_put_in_pot;
            players_stack[player_idx] -= amount_to_put_in_pot;
        } else if (action == Action::ALL_IN) {
            amount_to_put_in_pot = players_stack[player_idx];
            history_amount = amount_to_put_in_pot; // Amount added
            double total_contribution_this_action = pot_contribution_this_round[player_idx] + amount_to_put_in_pot;
            
            if (total_contribution_this_action > current_bet_to_match) { // All-in is a raise
                last_bet_or_raise_increment_this_round = total_contribution_this_action - current_bet_to_match;
                current_bet_to_match = total_contribution_this_action;
                last_aggressor_idx_this_round = player_idx;
            } // If all-in is a call, current_bet_to_match doesn't change by this player.
            
            pot_contribution_this_round[player_idx] = total_contribution_this_action;
            players_stack[player_idx] = 0;

        } else if (action == Action::BET_RAISE) {
            // chosen_amount_for_bet_raise is the TOTAL amount player wants their contribution to be this round
            if (chosen_amount_for_bet_raise <= current_bet_to_match || chosen_amount_for_bet_raise <= pot_contribution_this_round[player_idx]) {
                 // This indicates an invalid bet size passed to apply_action, should be caught by legal_actions
                 // Or it's an attempt to bet less than current commitment or call amount.
                 // For safety, convert to a CALL or CHECK if possible, or throw error.
                 // This should ideally not happen if legal_actions are followed.
                 // Let's assume chosen_amount_for_bet_raise is valid as per legal_actions.
            }

            amount_to_put_in_pot = chosen_amount_for_bet_raise - pot_contribution_this_round[player_idx];
            history_amount = amount_to_put_in_pot; // Amount added

            if (amount_to_put_in_pot >= players_stack[player_idx] - 1e-5) { // Bet/Raise is all-in
                amount_to_put_in_pot = players_stack[player_idx];
                history_amount = amount_to_put_in_pot;
                current_round_bets[player_idx] = Action::ALL_IN; // Treat as all-in
                action = Action::ALL_IN; // Update action for history
                chosen_amount_for_bet_raise = pot_contribution_this_round[player_idx] + players_stack[player_idx];
            }
            
            double old_bet_to_match_for_increment = current_bet_to_match; // Store current_bet_to_match before updating

            // Update state
            pot_contribution_this_round[player_idx] += amount_to_put_in_pot;
            players_stack[player_idx] -= amount_to_put_in_pot;
            
            // current_bet_to_match is the new level to be matched by others.
            current_bet_to_match = chosen_amount_for_bet_raise; 
            last_aggressor_idx_this_round = player_idx;

            // Update last_bet_or_raise_increment_this_round
            // This is the size of the bet or raise itself.
            // If it's the first bet in the round (old_bet_to_match_for_increment was effectively 0),
            // the increment is the total bet amount. Otherwise, it's the difference.
            if (old_bet_to_match_for_increment < 1e-5) { // Effectively 0
                last_bet_or_raise_increment_this_round = chosen_amount_for_bet_raise;
            } else { // It's a raise
                last_bet_or_raise_increment_this_round = chosen_amount_for_bet_raise - old_bet_to_match_for_increment;
            }
            // Ensure increment is not negative due to floating point issues
            if (last_bet_or_raise_increment_this_round < 0) last_bet_or_raise_increment_this_round = 0;
        }

        // Update cumulative pot after each action
        double current_round_pot = 0.0;
        for (double contribution : pot_contribution_this_round) {
            current_round_pot += contribution;
        }
        cumulative_pot = current_round_pot;
        for (double contribution : cumulative_pot_contribution) {
            cumulative_pot += contribution;
        }

        // Add to round_action_history *before* advancing round/player or checking for game over
        // Ensure action reflects if it became an ALL_IN
        Action effective_action = current_round_bets[player_idx]; // Use the potentially updated action (e.g. CALL -> ALL_IN)
        if (action != Action::DEAL && action != Action::POST_SB && action != Action::POST_BB) {
             round_action_history[round].emplace_back(player_idx, effective_action, history_amount);
        }

        if (game_over) { // e.g. player folded and game ends
             advance_round_if_needed_to_showdown();
             return;
        }
        
        advance_to_next_player(); // Tentatively advance, betting_round_complete will confirm.

        if (betting_round_complete()) {
            // Normal round advancement when betting round is complete
            advance_round();
            if (should_end_game()) {
                advance_round_if_needed_to_showdown();
                game_over = true;
            } else {
                // Next player for new round already set in advance_round()
                if(is_chance_node()) {
                    // The game loop should call apply_action(DEAL) next.
                }
            }
        } else {
            // Betting continues in the current round. next_player_idx is already set.
        }

        sync_attributes();
    }
    
    void advance_round_if_needed_to_showdown() {
        while(round != "showdown" && active_players_indices.size() > 1) {
            // Deal remaining community cards if not all players are all-in or folded
            bool all_remaining_are_all_in = true;
            if (active_players_indices.size() > 1) {
                for(int p_idx : active_players_indices) {
                    if (players_stack[p_idx] > 0) {
                        all_remaining_are_all_in = false;
                        break;
                    }
                }
            } else { // 0 or 1 active player, game ends.
                 all_remaining_are_all_in = true;
            }


            if (round == "preflop" && community_cards.size() < 3) deal_cards(); // Flop
            else if (round == "flop" && community_cards.size() < 4) deal_cards();   // Turn
            else if (round == "turn" && community_cards.size() < 5) deal_cards();    // River
            
            // If players are all-in, no more betting rounds, just advance state.
            if (all_remaining_are_all_in || active_players_indices.size() <=1 ) {
                 // Move contributions to cumulative
                for (int p = 0; p < game.num_players; p++) {
                    cumulative_pot_contribution[p] += pot_contribution_this_round[p];
                    pot_contribution_this_round[p] = 0.0;
                }
                current_bet_to_match = 0.0; // Reset for conceptual clarity, though no more betting
                last_bet_or_raise_increment_this_round = 0.0;
                last_aggressor_idx_this_round = -1;


                if (round == "preflop") round = "flop";
                else if (round == "flop") round = "turn";
                else if (round == "turn") round = "river";
                else if (round == "river") round = "showdown";
            } else { // Should not happen if game is ending.
                break; 
            }
        }
        game_over = true;
    }


    std::vector<double> returns() const {
        std::vector<double> rewards(game.num_players, 0.0);
        if (!game_over) return rewards;

        // Make a copy of the pot contributions to modify
        std::vector<double> cumulative_pot_contribution_copy = cumulative_pot_contribution;
        std::vector<double> pot_contribution_this_round_copy = pot_contribution_this_round;

        // Move all pot contributions to cumulative
        for (int p = 0; p < game.num_players; p++) {
            if (pot_contribution_this_round_copy[p] > 0) {
                 cumulative_pot_contribution_copy[p] += pot_contribution_this_round_copy[p];
                 pot_contribution_this_round_copy[p] = 0.0;
            }
        }

        double total_pot_value = 0.0;
        for (double x : cumulative_pot_contribution_copy) {
            total_pot_value += x;
        }
        
        // Rake not applied here, but if it were:
        // double rake_amount = total_pot_value * game.rake_percentage;
        // total_pot_value -= rake_amount;

        if (active_players_indices.size() == 1) {
            int winner_idx = *active_players_indices.begin();
            for (int p = 0; p < game.num_players; p++) {
                if (p == winner_idx) {
                    rewards[p] = total_pot_value - cumulative_pot_contribution_copy[p];
                } else {
                    rewards[p] = -cumulative_pot_contribution_copy[p];
                }
            }
        } else {
            // Multiple active players, need to evaluate hands
            PokerEvaluator evaluator;
            std::vector<int> hand_values;
            std::vector<int> active_players_vec(active_players_indices.begin(), active_players_indices.end());
            
            // Get hand values for active players
            for (int p : active_players_vec) {
                std::vector<Card> hole_cards = {cards[p*2], cards[p*2+1]};
                int hand_value = evaluator.evaluateHand(hole_cards, community_cards);
                hand_values.push_back(hand_value);
            }
            
            // Find best hand value
            int best_hand_value = *std::max_element(hand_values.begin(), hand_values.end());
            
            // Count winners (players with best hand value)
            int num_winners = std::count(hand_values.begin(), hand_values.end(), best_hand_value);
            
            // Split pot among winners
            double prize_per_winner = total_pot_value / num_winners;
            
            for (int p = 0; p < game.num_players; p++) {
                auto it = std::find(active_players_vec.begin(), active_players_vec.end(), p);
                if (it != active_players_vec.end()) {
                    int active_idx = std::distance(active_players_vec.begin(), it);
                    if (hand_values[active_idx] == best_hand_value) {
                        rewards[p] = prize_per_winner - cumulative_pot_contribution_copy[p];
                    } else {
                        rewards[p] = -cumulative_pot_contribution_copy[p];
                    }
                } else {
                    rewards[p] = -cumulative_pot_contribution_copy[p];
                }
            }
        }
        return rewards;
    }

    std::string to_string(bool show_cards = true) const {
        std::ostringstream oss;
        oss << "Round: " << round;
        oss << "\nCommunity Cards: ";
        for (const Card &c : community_cards) oss << c.toString() << " ";
        
        oss << "\nPlayers:\n";
        for (int p = 0; p < game.num_players; ++p) {
            oss << "  Player " << p << ": Stack=" << players_stack[p]
                << ", PotContRound=" << pot_contribution_this_round[p]
                << ", CumPotCont=" << cumulative_pot_contribution[p];
            if (cards.size() >= static_cast<size_t>((p + 1) * 2)) {
                 oss << ", Cards=" << cards[p*2].toString() << cards[p*2+1].toString();
            }
            if (active_players_indices.count(p)) oss << " (Active)";
            oss << "\n";
        }
        oss << "Current Bet to Match: " << current_bet_to_match << "\n";
        oss << "Last Raise Inc: " << last_bet_or_raise_increment_this_round << "\n";
        oss << "Last Aggressor: " << last_aggressor_idx_this_round << "\n";
        oss << "Next Player: " << next_player_idx << "\n";
        oss << "Status: " << (game_over ? "Game Over" : "In Progress") << "\n";
        oss << "Is Chance Node: " << (is_chance_node() ? "Yes" : "No") << "\n";
        oss << "\nRound Action History:\n";
        for (const auto& round_entry : round_action_history) {
            oss << "  " << betting_round_to_string(round_entry.first) << ": ";
            for (const auto& action_tuple : round_entry.second) {
                oss << "[P" << std::get<0>(action_tuple)
                    << " " << action_to_string(std::get<1>(action_tuple), std::get<2>(action_tuple)) // Pass amount here
                    << "] ";
            }
            oss << "\n";
        }
        return oss.str();
    }
    
    std::string betting_round_to_string(const std::string& br) const {
        return br;  // Just return the string directly since it's already in string format
    }

    std::string get_player_hole_cards_str(int player_idx, bool full_info = true) const {
        if (!full_info && player_idx != current_player() && !game_over) return "??"; // Hide opponent cards unless showdown/current player
        if (cards.size() >= static_cast<size_t>((player_idx + 1) * 2)) {
            return cards[player_idx * 2].toString() + cards[player_idx * 2 + 1].toString();
        }
        return "";
    }

    std::string get_community_cards_str() const {
        std::string s;
        for (const auto& card : community_cards) {
            s += card.toString() + " ";
        }
        return s;
    }

    // Helper to get total pot size
    double get_total_pot() const {
        double total = 0;
        for(double c : cumulative_pot_contribution) total += c;
        for(double c : pot_contribution_this_round) total += c;
        return total;
    }

    // Helper method to calculate SPR after a bet is called
    double calculate_spr_after_call(double bet_total_contribution, int opponent_idx) const {
        double pot_after_call = get_total_pot() + (bet_total_contribution - pot_contribution_this_round[current_player()]);
        return players_stack[opponent_idx] / pot_after_call;
    }

    // Helper method to merge similar bet sizes based on merging threshold
    std::vector<double> merge_bet_sizes(std::vector<double> bet_sizes, double merging_threshold) const {
        if (bet_sizes.empty()) return bet_sizes;
        
        // Sort in descending order
        std::sort(bet_sizes.rbegin(), bet_sizes.rend());
        
        std::vector<double> merged_sizes;
        std::vector<bool> removed(bet_sizes.size(), false);
        
        double current_pot = get_total_pot();
        
        for (size_t i = 0; i < bet_sizes.size(); ++i) {
            if (removed[i]) continue;
            
            double x_percent = (bet_sizes[i] / current_pot) * 100.0;
            merged_sizes.push_back(bet_sizes[i]);
            
            // Remove similar sized bets
            for (size_t j = i + 1; j < bet_sizes.size(); ++j) {
                if (removed[j]) continue;
                
                double y_percent = (bet_sizes[j] / current_pot) * 100.0;
                if ((100.0 + x_percent) / (100.0 + y_percent) < 1.0 + merging_threshold) {
                    removed[j] = true;
                }
            }
        }
        
        return merged_sizes;
    }
};


// Operator overload for printing legal actions (optional, for debugging)
std::ostream& operator<<(std::ostream& os, const std::vector<std::pair<Action, double>>& actions) {
    os << "[";
    for (size_t i = 0; i < actions.size(); ++i) {
        os << action_to_string(actions[i].first, actions[i].second);
        if (i < actions.size() - 1) os << ", ";
    }
    os << "]";
    return os;
}
std::ostream& operator<<(std::ostream& os, const std::vector<double>& rewards) {
    os << "[";
    for (size_t i = 0; i < rewards.size(); ++i) {
        os << rewards[i];
        if (i < rewards.size() - 1) os << ", ";
    }
    os << "]";
    return os;
}


// // Example main for testing (optional)
// int main() {
//     HUGame game;
//     HUState state = game.new_initial_state();
    
//     std::cout << "Starting new Heads-Up Poker game!\n";
//     std::cout << "Initial State:\n" << state.to_string(true) << std::endl;

//     int max_actions = 100;
//     int actions_taken = 0;

//     while (!state.game_over && actions_taken < max_actions) {
//         actions_taken++;
//         std::cout << "\n-----------------------------------\n";
//         std::cout << "Round: " << state.round << "\n";
//         std::cout << "Community Cards: " << state.get_community_cards_str() << "\n";
        
//         // Only show current player's cards unless game is over
//         for (int p = 0; p < game.num_players; p++) {
//             std::cout << "Player " << p << " Cards: " << state.get_player_hole_cards_str(p, false) << "\n";
//         }
        
//         std::cout << "Current Bet: " << state.current_bet << "\n";
//         std::cout << "Player 0 Stack: " << state.players_stack[0] << "\n";
//         std::cout << "Player 1 Stack: " << state.players_stack[1] << "\n";
//         std::cout << "Total Pot: " << state.get_total_pot() << "\n";

//         std::vector<std::pair<Action, double>> legal = state.legal_actions(game.pot_fraction_bet_sizes, game.fixed_bet_sizes_bb);

//         if (state.is_chance_node()) {
//             std::cout << "Dealing cards...\n";
//             state.apply_action(Action::DEAL);
//             continue;
//         }

//         std::cout << "\nPlayer " << state.current_player() << "'s turn. Legal actions:\n";
//         for (size_t i = 0; i < legal.size(); i++) {
//             std::cout << i << ": " << action_to_string(legal[i].first, legal[i].second) << "\n";
//         }

//         int choice;
//         do {
//             std::cout << "Enter the number of your chosen action: ";
//             std::cin >> choice;
//         } while (choice < 0 || choice >= legal.size());

//         Action chosen_action = legal[choice].first;
//         double chosen_amount = legal[choice].second;

//         try {
//             if (chosen_action == Action::BET_RAISE) {
//                 state.apply_action(chosen_action, chosen_amount);
//             } else {
//                 state.apply_action(chosen_action);
//             }
//         } catch (const std::runtime_error& e) {
//             std::cerr << "Error applying action: " << e.what() << std::endl;
//             break;
//         }
//     }

//     std::cout << "\n===================================\n";
//     std::cout << "Game Over!\n";
//     std::cout << "Final State:\n" << state.to_string(true) << std::endl;
//     std::cout << "Returns: " << state.returns() << std::endl;

//     return 0;
// }


// After HUState is fully defined, implement HUGame::new_initial_state
HUState HUGame::new_initial_state() const {
     return HUState(*this);
}

HUGame::HUGame(
    int num_players_val,
    double initial_stack_val,
    double small_blind_amount_val,
    double big_blind_amount_val,
    float rake_percentage_val,
    const std::vector<double>& pot_fractions,
    const std::vector<double>& fixed_bbs,
    double all_in_threshold_val,
    double force_all_in_threshold_val,
    double merging_threshold_val
) : num_players(num_players_val),
    initial_stack(initial_stack_val),
    small_blind_amount(small_blind_amount_val),
    big_blind_amount(big_blind_amount_val),
    rake_percentage(rake_percentage_val),
    pot_fraction_bet_sizes(pot_fractions),
    fixed_bet_sizes_bb(fixed_bbs),
    all_in_threshold(all_in_threshold_val),
    force_all_in_threshold(force_all_in_threshold_val),
    merging_threshold(merging_threshold_val)
{
    // Constructor body can be empty as we're using initializer list
}
