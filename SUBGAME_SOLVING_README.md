# Subgame MCCFR Implementation for Heads-Up No-Limit Texas Hold'em

## Overview

This implementation extends the existing Monte-Carlo CFR solver with **Subgame CFR** capability, allowing you to solve specific situations in Heads-Up No-Limit Texas Hold'em starting from custom states rather than the full game tree.

## Key Features

### 1. **Subgame Initialization**
- Start from any predefined game state
- Define betting action history (e.g., "c-r200-c" for check-raise200-call)
- Set community cards for any street (flop, turn, or river)

### 2. **Chance Node Handling**
- Unknown hole cards for both players are sampled uniformly
- Remaining community cards are sampled from the deck
- Proper integration with CFR computation

### 3. **Game Tree Construction**
- Only explores forward from the specified state
- Simulates all legal actions for current and future streets
- Supports all bet sizing options available in the base game

## Usage

### Interactive Mode

Run the solver and choose option for "Solve subgame from custom state (Subgame CFR)":

```bash
./solver_cfr_rt
```

You'll be prompted to enter:
1. **Betting History**: Format like "c-r200-c" (check-raise200-call) or empty for start of round
2. **Community Cards**: Format like "As Kd Qh" (A♠ K♦ Q♥) or empty for preflop
3. **Number of CFR iterations**: How many iterations to run
4. **Number of threads**: For parallel processing

### Programmatic Usage

```cpp
#include "solver_cfr_rt.cpp"

// Initialize game and clusters
HUGame game_config;
preloadClusters_ic();

// Define subgame parameters
std::string betting_history = "c-r200-c";        // Check-raise200-call
std::string community_cards = "As Kd Qh";        // Flop: A♠ K♦ Q♥  
int iterations = 1000;
std::string output_file = "subgame_strategy.csv";
int num_threads = 4;

// Solve the subgame
solve_subgame_cfr(game_config, betting_history, community_cards, 
                  iterations, output_file, num_threads);
```

## Input Formats

### Betting History
- Actions separated by dashes: `"-"`
- **Check**: `"c"` or `"check"`
- **Call**: `"c"`
- **Fold**: `"f"` or `"fold"`
- **Bet**: `"b100"` (bet 100 chips)
- **Raise**: `"r200"` (raise to 200 chips)
- **Examples**: 
  - `"c-r200-c"` = check, raise to 200, call
  - `"b100-c"` = bet 100, call
  - `""` = no betting history (start of round)

### Community Cards
- Cards separated by spaces
- Format: `[Rank][Suit]` where:
  - **Ranks**: `2, 3, 4, 5, 6, 7, 8, 9, 10, J, Q, K, A`
  - **Suits**: `c, d, h, s` (clubs, diamonds, hearts, spades)
- **Examples**:
  - `"As Kd Qh"` = A♠ K♦ Q♥ (flop)
  - `"As Kd Qh 7c"` = A♠ K♦ Q♥ 7♣ (turn)
  - `"As Kd Qh 7c 2s"` = A♠ K♦ Q♥ 7♣ 2♠ (river)
  - `""` = no community cards (preflop)

## Algorithm Details

### Subgame CFR Process

1. **State Initialization**: Parse betting history and community cards to create the starting state
2. **Hole Card Sampling**: For each CFR iteration, randomly sample hole card combinations for both players
3. **Chance Node Resolution**: Sample remaining community cards when needed (e.g., turn/river cards)
4. **CFR Computation**: Apply standard CFR algorithm from the subgame root
5. **Strategy Aggregation**: Accumulate regrets and strategies across all sampled scenarios

### Key Differences from Full Game CFR

- **Reduced Game Tree**: Only explores forward from the specified state
- **Sampling Integration**: Handles unknown private information through sampling
- **Information Set Abstraction**: Uses the same clustering system as full game CFR
- **Memory Efficiency**: Smaller game tree means lower memory requirements

## Example Scenarios

### Scenario 1: Flop Decision Point
```
Betting History: "c-r200-c"     # Check-raise200-call
Community Cards: "As Kd Qh"     # A♠ K♦ Q♥
```
Solves optimal play from this flop situation where there's been action.

### Scenario 2: Turn All-In Spot  
```
Betting History: "c-c-b150"     # Check-check-bet150
Community Cards: "As Kd Qh 7c"  # A♠ K♦ Q♥ 7♣
```
Analyzes turn play after flop was checked through and turn was bet.

### Scenario 3: River Showdown
```
Betting History: ""              # Start of betting round
Community Cards: "As Kd Qh 7c 2s" # A♠ K♦ Q♥ 7♣ 2♠
```
Solves river play with all community cards revealed.

## Files Generated

- **Strategy CSV**: Contains optimal strategies for each information set
  - Format: `infoset_key, action_string, probability`
  - Information sets include: round, player, pot size, hand abstraction
  - Actions include: check, call, fold, bet/raise with amounts

## Performance Notes

- **Iterations**: Start with 1000-10000 for testing, scale up for production
- **Threads**: Use up to your CPU core count for optimal performance  
- **Memory**: Subgame CFR uses significantly less memory than full game CFR
- **Convergence**: Smaller subgames typically converge faster than full games

## Integration with Existing Code

The subgame CFR implementation extends the existing CFR solver without breaking compatibility:
- Uses same clustering system (`preloadClusters_ic()`)
- Same information set generation (`get_infoset_key_for_state()`)
- Same action representation and game logic
- Compatible with existing strategy analysis tools

## Testing

Run the test program to verify functionality:

```bash
g++ -std=c++17 -pthread -O3 test_subgame_cfr.cpp -o test_subgame_cfr
./test_subgame_cfr
```

This will run three test scenarios and generate strategy files for verification. 
