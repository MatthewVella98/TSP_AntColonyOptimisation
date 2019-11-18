# Travelling Salesman Problem 

### Abstract

TSP asks the following question: "Given a list of cities and the distances between each pair of cities, what is the shortest possible route that visits each city and returns to the origin city?". This is an Np-hard problem. 


The problem was first studied in 1930 and is one of the most intensively studied problems in optimisation. This problem is used as a benchmark for many optimisation methods.  Many algorithms exist (heuristics and exact) which solve this problem (Exact algorithms are obviously not optimal).


The TSP has several applications: Planning logistics, DNA, Astronomy etc. 


### Description 

#### TSP as a Graph problem

The TSP, can be modelled as un undirected weighted graph, such that:
  - The cities are the graph's vertices
  - The paths are the graph's edges
  - A path's istance is the edge's weight 

This is a minimisation problem, which starts and finishes at a particular vertes, after visiting each vertex exactly once. The model is very often a complete graph. 

### Asymmetric and Symmetric 

Symmetric TSP: The distance between two cities is the same in each opposite direction, forming an undirected graph. This symmetry halves the number of possible solutions. 


Asymmetric TSP: Paths may not exist in both directions/, forming a directed graph.
