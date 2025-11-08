# Pathfinding Visualization with Raylib and C++

## Project Overview

This is a educational C++ project that implements a pathfinding algorithm visualization tool using Raylib. The application allows users to interactively create a grid-based environment and demonstrate Dijkstra's pathfinding algorithm.

## Features

- **Interactive Grid Visualization**: Create a grid of squares using Raylib
- **Node Selection**:
  - <b>Green Key</b>: Mark source node
  - <b>Red Key</b>: Mark destination node
  - <b>M Key</b>: Create obstacles (magenta cells)
- **Pathfinding**:
  - <b>Shift+G</b>: Find the shortest path using Dijkstra's algorithm
  - Searched spaces are marked in dark green
  - Shortest path is highlighted in blue
- **Exit**:
  - <b>Escape</b> or <b>Q Key</b>: Quit the application

## Learning Objectives

- Practice C++ programming
- Implement Dijkstra's pathfinding algorithm
- Gain experience with Raylib graphics library

## Planned Improvements

- Implement A* pathfinding algorithm
- Add more advanced obstacle creation methods
- Enhance visualization features

## Prerequisites

- C++ compiler
- Raylib library
- Basic understanding of pathfinding algorithms

## Compilation

Provide specific compilation instructions for your project (e.g., CMake, Makefile, or specific compiler flags)

## Usage

1. Launch the application
2. Use keys to interact with the grid:
   - Green Key: Set source node
   - Red Key: Set destination node
   - M Key: Create obstacles
   - Shift+G: Find path
   - Escape/Q: Quit

## Algorithm Details

Currently implements <b>Dijkstra's algorithm</b> for finding the shortest path between two points on a grid, with obstacles.

## Future Work

- Implement A* pathfinding algorithm
- Add more advanced visualization techniques
- Improve performance and add more interactive features
