# Texas Hold'em Abstraction with Rake

A comprehensive implementation of Heads-Up No-Limit Texas Hold'em with card abstraction, Monte Carlo Counterfactual Regret Minimization (MCCFR), and subgame solving capabilities.

## Overview

This project implements a complete Texas Hold'em poker engine with advanced solving techniques. It features:

- **Full Game Tree Solver**: MCCFR implementation for complete game analysis
- **Subgame Solver**: Solve specific game situations from custom states
- **Card Abstraction**: Hand clustering for efficient computation
- **Multi-threaded Processing**: Parallel computation for faster convergence
- **Exploitability Tracking**: Monitor solution quality during training
- **Real-time Strategy Access**: Query optimal strategies during gameplay

## Features

### Core Game Engine
- **Heads-Up No-Limit Texas Hold'em** with proper betting rules
- **Blind posting** and **position-based** play order
- **All betting actions**: check, call, fold, bet, raise, all-in
- **Multiple bet sizing options**: pot fractions, fixed BB amounts
- **Proper hand evaluation** with all poker hand rankings
- **Rake calculation** support (configurable percentage)

### AI Solving Capabilities
- **Monte Carlo CFR**: Efficient regret minimization with sampling
- **Outcome Sampling**: Accelerated convergence for large games
- **Linear CFR**: Advanced regret minimization techniques
- **Subgame Solving**: Analyze specific game situations
- **Exploitability Computation**: Measure solution quality
- **Strategy Export/Import**: Save and load computed strategies

### Card Abstraction System
- **Hand Clustering**: Group similar hands for computational efficiency
- **Multi-street Support**: Preflop, flop, turn, and river abstractions
- **Equity-based Clustering**: Hands grouped by similar expected values
- **Configurable Granularity**: Balance between accuracy and speed

## Installation

### Prerequisites
- **C++17** compatible compiler (GCC 7+, Clang 5+, MSVC 2017+)
- **Make** or **CMake** for building
- **Python 3.6+** (for utility scripts, optional)

### Building from Source

1. **Clone the repository**:
   ```bash
   git clone <repository-url>
   cd texas_holdem_abstraction_with_rake
   ```

2. **Download Required CSV Files**:
   
   **⚠️ IMPORTANT**: You must download the necessary CSV files for the card abstraction system to work properly. These files contain pre-computed equity data and cluster mappings.
   
   Download from: [https://www.dropbox.com/scl/fo/ln3e95uryi8skr46b1v7y/ABgAZgyl1Pu0VJpHhVOqSr8?rlkey=7ajiuhns3vxxomliyhq03pqh7&st=e1bb272d&dl=0](https://www.dropbox.com/scl/fo/ln3e95uryi8skr46b1v7y/ABgAZgyl1Pu0VJpHhVOqSr8?rlkey=7ajiuhns3vxxomliyhq03pqh7&st=e1bb272d&dl=0)
   
   Place the downloaded files in the `utils/` folder:
   ```
   utils/
   ├── repr_l2_flop_equities_clustered_hands_with_avg_equity.csv
   ├── repr_l2_turn_equities_clustered_hands_with_avg_equity.csv
   └── repr_l2_river_equities_clustered_hands_with_avg_equity.csv
   ```

3. **Build the project**:
   ```bash
   # Using Make
   make
   
   # Or manually with g++
   g++ -std=c++17 -pthread -O3 realtime_solver.cpp -o solver_cfr_rt
   ```

## Usage

### Basic Game Solving

Solve a complete heads-up game from the beginning:

```bash
./solver_cfr_rt
```

Choose option 1: "Solve full game from initial state (MCCFR)"

### Subgame Solving

Solve specific game situations:

```bash
./solver_cfr_rt
```

Choose option 2: "Solve subgame from custom state (Subgame CFR)"

You'll be prompted for:
- **Betting History**: e.g., "c-r200-c" (check-raise200-call)
- **Community Cards**: e.g., "As Kd Qh" (A♠ K♦ Q♥)
- **Iterations**: Number of CFR iterations
- **Threads**: Number of parallel threads

### Programmatic Usage

```cpp
#include "hu.cpp"
#include "realtime_solver.cpp"

// Initialize game configuration
HUGame game_config(2, 100.0, 0.5, 1.0, 0.07);

// Load card abstractions
preloadClusters_ic();

// Solve full game
solve_from_state_mccfr_accelerated(
    game_config.new_initial_state(),
    10000,                    // iterations
    "strategy.csv",           // output file
    4,                        // threads
    true,                     // use outcome sampling
    true,                     // use linear CFR
    100                       // exploitability check interval
);

// Solve subgame
solve_subgame_cfr(
    game_config,
    "c-r200-c",              // betting history
    "As Kd Qh",              // community cards
    5000,                    // iterations
    "subgame_strategy.csv",  // output file
    4                        // threads
);
```

## Configuration

### Game Parameters

```cpp
HUGame game(
    2,                    // num_players
    100.0,               // initial_stack (BB)
    0.5,                 // small_blind_amount
    1.0,                 // big_blind_amount
    0.07,                // rake_percentage (7%)
    {0.33, 0.5, 0.75, 1.0},  // pot_fraction_bet_sizes
    {2.5, 3.0, 4.0},         // fixed_bet_sizes_bb
    1.5,                 // all_in_threshold
    0.2,                 // force_all_in_threshold
    0.1                  // merging_threshold
);
```

### Betting Options

- **Pot Fraction Bets**: 33%, 50%, 75%, 100% of pot
- **Fixed BB Bets**: 2.5BB, 3BB, 4BB
- **All-in Thresholds**: Configurable for different game dynamics

## File Structure

```
├── hu.cpp                    # Core game engine (HUState, HUGame, PokerEvaluator)
├── realtime_solver.cpp       # MCCFR solver and utilities
├── hu_demo.cpp              # Demo program
├── utils/                   # Card abstraction CSV files
│   ├── repr_l2_flop_equities_clustered_hands_with_avg_equity.csv
│   ├── repr_l2_turn_equities_clustered_hands_with_avg_equity.csv
│   └── repr_l2_river_equities_clustered_hands_with_avg_equity.csv
├── build/                   # Build artifacts
├── SUBGAME_SOLVING_README.md # Detailed subgame solving documentation
└── README.md               # This file
```

## Key Components

### HUState Class
- **Game State Management**: Tracks current round, pot, player stacks, actions
- **Action Validation**: Ensures all actions are legal according to poker rules
- **Round Progression**: Handles preflop → flop → turn → river → showdown
- **Pot Management**: Tracks contributions, side pots, and all-in scenarios

### PokerEvaluator Class
- **Hand Ranking**: Evaluates 5-card poker hands
- **Combination Generation**: Finds best 5-card combination from 7 cards
- **Score Calculation**: Numeric hand strength for comparison

### MCCFR Solver
- **Regret Minimization**: Implements counterfactual regret minimization
- **Information Set Abstraction**: Groups similar game states
- **Strategy Computation**: Calculates optimal action probabilities
- **Convergence Tracking**: Monitors solution quality over iterations

## Performance Considerations

### Memory Usage
- **Card Abstractions**: ~1GB for all CSV files
- **Strategy Storage**: Scales with number of information sets
- **Game Tree**: Significantly reduced with abstractions

### Computation Time
- **Full Game**: 10K-100K iterations for convergence
- **Subgames**: 1K-10K iterations typically sufficient
- **Parallelization**: Linear speedup with thread count

### Optimization Tips
- Use **outcome sampling** for large games
- Enable **linear CFR** for faster convergence
- **Tune abstraction granularity** based on accuracy needs
- **Monitor exploitability** to determine convergence

## Examples

### Example 1: Basic Game Solving
```bash
./solver_cfr_rt
# Choose: 1 (Full game)
# Iterations: 10000
# Threads: 4
# Output: strategy.csv
```

### Example 2: Flop Decision Analysis
```bash
./solver_cfr_rt
# Choose: 2 (Subgame)
# Betting: c-r200-c
# Community: As Kd Qh
# Iterations: 5000
# Threads: 4
# Output: flop_analysis.csv
```

### Example 3: River All-in Spot
```bash
./solver_cfr_rt
# Choose: 2 (Subgame)
# Betting: b150-a
# Community: As Kd Qh 7c 2s
# Iterations: 3000
# Threads: 4
# Output: river_allin.csv
```

## Troubleshooting

### Common Issues

1. **CSV Files Missing**: Ensure all required CSV files are in `utils/` folder
2. **Compilation Errors**: Verify C++17 support and required libraries
3. **Memory Issues**: Reduce thread count or use smaller abstractions
4. **Slow Convergence**: Increase iterations or enable outcome sampling

### Debug Mode

Compile with debug flags for detailed output:
```bash
g++ -std=c++17 -pthread -O0 -g -DDEBUG realtime_solver.cpp -o solver_cfr_rt_debug
```

## Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests if applicable
5. Submit a pull request

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Citation

If you use this code in research, please cite:
```bibtex
@software{aof_mccfr,
  title={HeadsUP Texas-Holdem real-time MCCFR Solver: High-Performance Game Theory Optimal Poker Training},
  author={Zakaria El Jaafari},
  year={2025},
  url={https://github.com/nier2kirito/AoF_mccfr}
}
```
## Acknowledgments

- **CFR Algorithm**: Based on the counterfactual regret minimization framework
- **Card Abstractions**: Pre-computed equity data for efficient computation
- **Poker Rules**: Implementation follows standard Texas Hold'em rules

## Support

For questions and issues:
- Check the existing documentation
- Review the code comments
- Open an issue on the repository
- Contact the maintainers

---

**Note**: This system requires the CSV files containing clustered hand representations for proper operation. Make sure to download them from the provided link and place them in the `utils/` folder before running the solver.
