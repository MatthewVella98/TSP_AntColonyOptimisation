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

### Integer linear programming formulations 

  - Miller-Tucker-Zemlin formulation
  - Dantzig-Fulkerson-Johnson formulation

### Computing a solution 

NP-hard problms are usually solved in the following ways:
  - Building exact algorithms: These work fast only for small problem sizes. 
  - Building suboptimal or heuristic algorithms: These algorithms deliver good qualiy            approximated solutions in a reasonable time.
  - Finding sub-problems or special cases for which either better or exact heuristics are        possible.

#### Exact algorithms 

  - Brute-force search: Trying out all permutations and checking out which one is the cheapest. The problem with this search method is that the running time is o(n!), the factorial of the number of cities. With 5 cities, the running time is reasoanble o(5!) = 120 searches. But if the number of cities is pumped up to 15, the running time is 0(15!) = 1.3076744e+12 searches, for just 15 cities. 

   - Dynamic programming approach - Held-Karp algorithm: This algorithm computes the solutions of all subproblem starting with the smallest. Whenever computing a solution requires solutions for smaller, look up the solutions which are already computed. This is solved recursively. 

Other approaches include:
   - Branch-and-bound algorithms, which are used to process TSPs which contain 40-60 cities. 
   - Progressive improvement algorithms: These use techniques similar to linear programming      and works well for up to 200 cities.
   - Implementations of branch-and-bound and branch-an-cut: This method should be used to         solve a very large number of instances. Currently, this method holds the record with      85,900 cities.  