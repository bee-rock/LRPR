LRPR
====

CPSC 517: Sparse Matrix Computation
Low Rank Page Rank: A course project in sparse matrix computation
Brock Hargreaves

The problem of Pagerank is a simple one to state: Given a collection of websites, how do we
rank them? The primary way of formulating this utilizes a transition matrix which relates how web pages interact with each other.

We investigate what the effect of a low rank approximation for the transition matrix has on the power method and an inner-outer iteration for solving the Pagerank problem.

The purpose of the low rank approximation is two fold: (1) to reduce memory requirements (2) to decrease computational time. We show that we see an improvement in storage requirements and a decrease in computational time if we discard the time it takes to perform the low rank approximation, however at the sacrifice of accuracy.

-----------------------------

Folders:

         doc: PDF's for project and presentation and associated papers
   utilities: Various tools use for solving the pagerank problem
    examples: Examples from project and presentation
innout-small: Code and data from: "An inner-outer iteration for computing 
              pagerank"

Add the utilities folder to your MATLAB path and execute 
project_examples.m. Note that you will need to  make a small edit to 
the working directory string and example string if you want.


Notes:

1. The utilities priorityqueue and bigraph may need to be recompiled to 
be compatible with your system. Simply run their demo/test files.

2. I wrote my own implementations for each of the algorithms. However,
without a lot optimization in mind. One should try using the algorithms
included in innout-small.
