#include <iostream>
#include <vector>
#include <string>
#include <random>
#include <stdexcept>
#include <algorithm> // For std::min, std::find
#include <iomanip>   // For std::fixed, std::setprecision if needed for printing

// Include the entire hu.cpp. For a larger project, you'd typically use header files.
// Ensure hu.cpp does not have its own main() function uncommented.
#include "hu.cpp" 

// Forward declaration if operator<< for actions is not picked up via include,
// but it should be if hu.cpp is included directly.
// std::ostream& operator<<(std::ostream& os, const std::vector<Action>& actions);

int main() {
    HUGame game;
    HUState state = game.new_initial_state();
    std::mt19937 rng(std::random_device{}()); // Random number generator

    std::cout << std::fixed << std::setprecision(2); // For consistent double printing

    std::cout << "Initial State:\n" << state.to_string() << std::endl;

    int max_actions_per_game = 50; // Safety break for the simulation loop
    int actions_taken_count = 0;

    while (!state.game_over && actions_taken_count < max_actions_per_game) {
        actions_taken_count++;
        std::cout << "\n-----------------------------------\n";
        std::cout << "Game Turn: " << actions_taken_count << "\n";
        std::cout << state.to_string();
        
        // Write infoset data for the current player
        // state.write_infoset_data(std::cout); 
        
        std::vector<Action> legal_actions_list = state.legal_actions();
        std::cout << "Legal actions for player " << state.current_player() << ": " << legal_actions_list << std::endl;

        if (state.is_chance_node()) {
            std::cout << "Player " << state.current_player() << " (Chance Node) chooses DEAL" << std::endl;
            state.apply_action(Action::DEAL);
            continue;
        }
        
        if (legal_actions_list.empty()) {
            std::cout << "No legal actions available for player " << state.current_player() << ". Game might be stuck or over." << std::endl;
            // This case should ideally be prevented by game logic marking game_over.
            if (!state.game_over) { // If game isn't marked over, force it to evaluate returns
                state.game_over = true; 
                state.advance_round_if_needed_to_showdown(); // Ensure it goes to showdown if appropriate
            }
            break;
        }

        // Simple random agent logic
        std::uniform_int_distribution<int> dist(0, legal_actions_list.size() - 1);
        Action chosen_action_type = legal_actions_list[dist(rng)];
        
        Action action_to_apply = chosen_action_type;
        double total_contribution_for_bet_raise = 0.0;

        std::cout << "Player " << state.current_player() << " considers " << action_to_string(chosen_action_type);

        if (chosen_action_type == Action::BET_RAISE) {
            double min_total_contrib_target = state.get_min_bet_raise_total_amount();
            double max_total_contrib_target = state.get_max_bet_raise_total_amount(); // Player's stack + current round pot contribution

            if (min_total_contrib_target >= max_total_contrib_target && max_total_contrib_target > state.pot_contribution_this_round[state.current_player()]) { 
                // Min valid bet/raise requires going all-in (and player has chips to bet).
                action_to_apply = Action::ALL_IN;
                std::cout << " (effectively ALL_IN as min_bet_target=" << min_total_contrib_target 
                          << " >= max_bet_target=" << max_total_contrib_target << ")";
            } else if (min_total_contrib_target < max_total_contrib_target) {
                // Can make a bet/raise that isn't necessarily all-in.
                // Agent policy: try to bet/raise to min_total_contrib_target, or min + BB, capped by all-in.
                total_contribution_for_bet_raise = min_total_contrib_target;
                
                // Example of a slightly more dynamic bet:
                if (min_total_contrib_target + BIG_BLIND_AMOUNT <= max_total_contrib_target) {
                    total_contribution_for_bet_raise = min_total_contrib_target + BIG_BLIND_AMOUNT;
                }
                // Ensure it doesn't exceed max (all-in amount)
                total_contribution_for_bet_raise = std::min(total_contribution_for_bet_raise, max_total_contrib_target);
                // Ensure it's at least the minimum valid bet
                total_contribution_for_bet_raise = std::max(total_contribution_for_bet_raise, min_total_contrib_target);


                std::cout << " with total target contribution: " << total_contribution_for_bet_raise;
                action_to_apply = Action::BET_RAISE; // Confirming it's a BET_RAISE
            } else {
                 // This case implies min_total_contrib_target >= max_total_contrib_target but player has no chips to bet (max_total_contrib_target == current pot contrib)
                 // Or min_total_contrib_target is somehow invalid.
                 // Fallback: if BET_RAISE was legal but logic here is tricky, try ALL_IN if stack > 0, or CHECK if available.
                if (state.players_stack[state.current_player()] > 0 &&
                    std::find(legal_actions_list.begin(), legal_actions_list.end(), Action::ALL_IN) != legal_actions_list.end()) {
                    action_to_apply = Action::ALL_IN;
                    std::cout << " (fallback to ALL_IN)";
                } else if (std::find(legal_actions_list.begin(), legal_actions_list.end(), Action::CHECK) != legal_actions_list.end()) {
                    action_to_apply = Action::CHECK;
                    std::cout << " (fallback to CHECK)";
                } else {
                    // Should not happen if legal_actions is not empty. Take first available.
                    action_to_apply = legal_actions_list[0]; 
                    std::cout << " (fallback to first legal: " << action_to_string(action_to_apply) << ")";
                }
            }
        }
        std::cout << " -> Applying: " << action_to_string(action_to_apply) << std::endl;
        
        try {
            if (action_to_apply == Action::BET_RAISE) {
                state.apply_action(action_to_apply, total_contribution_for_bet_raise);
            } else { 
                // Covers ALL_IN, FOLD, CHECK, CALL, DEAL (though DEAL handled above)
                state.apply_action(action_to_apply);
            }
        } catch (const std::runtime_error& e) {
            std::cerr << "Error applying action: " << e.what() << std::endl;
            std::cerr << "Current state before error:\n" << state.to_string() << std::endl;
            std::cerr << "Chosen action type: " << action_to_string(chosen_action_type) << ", Applied action: " << action_to_string(action_to_apply) << std::endl;
            if(action_to_apply == Action::BET_RAISE) std::cerr << "Bet amount: " << total_contribution_for_bet_raise << std::endl;
            break; // End game on error
        }
    }

    if (actions_taken_count >= max_actions_per_game && !state.game_over) {
        std::cout << "\nMax actions reached, game might be in a loop. Forcing end." << std::endl;
        state.game_over = true; // Force game over
        state.advance_round_if_needed_to_showdown(); // Try to resolve to showdown
    }

    std::cout << "\n===================================\n";
    std::cout << "Final State:\n" << state.to_string() << std::endl;
    std::cout << "Returns: " << state.returns() << std::endl;

    return 0;
}
