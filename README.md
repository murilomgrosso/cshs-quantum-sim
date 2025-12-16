# CSHS Quantum Simulation
A simple C++ quantum simulation of the CHSH game and its classical and quantum strategies.

This project demonstrates how quantum strategies outperform classical ones in a well-known theoretical game from quantum information theory.

## The CSHS Game
Theoretical framework used to experimentally test concepts from quantum information.

The game involves three entities:
- Alice (player)
- Bob (player)
- Referee

Rules:
1. The referee randomly chooses two bits:
    - `x` for Alice
    - `y` for Bob
2. Alice and Bob receive their respective bits without communicating.
3. Alice outputs a bit `a`, and Bob outputs a bit `b`.
4. Alice and Bob win if:
```
a ⊕ b = x ∧ y
```
Otherwise, they lose.

## Classical vs Quantum Strategies
- Classical strategies have a maximum winning probability of 75%.
- Quantum strategies, using entangled states, can achieve a winning probability of 85%.

This simulation compares both approaches and highlights the advantage provided by quantum mechanics.

```
1) Random Strategy 
2) Copy Strategy 
3) Opposite Strategy 
4) Always zero Strategy 
5) Always one Strategy 
6) Quantum Strategy 

WIN RATE (134000 games played)
1) [##########################........................] 50.1507%
2) [#############.....................................] 24.9358%
3) [#############.....................................] 24.9358%
4) [######################################............] 75.0194%
5) [######################################............] 75.0194%
6) [###########################################.......] 85.3398%
```

## How to build and run
> g++ src/*.cpp -I includes -o bin/main && bin/main

## References
[1] Riccardo Bassoli, Holger Boche, Christian Deppe, et al.: Quantum Communication Networks, Springer, 2021.
