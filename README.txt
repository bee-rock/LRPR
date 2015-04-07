LRPR
====

CPSC 517: Sparse Matrix Computation
Low Rank Page Rank: A course project in sparse matrix computation
Brock Hargreaves

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
