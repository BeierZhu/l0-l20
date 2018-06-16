# l0-l20 Toolbox
l_0-minimization or l_20-minimization solver

## Usage

### OMP.m

Orthogonal Matching Pursuit (OMP) solver for l_0 minimization problem:

minimize     ||x||_1
subject to   y = Ax,

function [x, Out] = OMP(A,y,varargin)

*Required inputs*

A -- an m x n matrix

y -- an m x 1 vector

*Optional inputs*

'maxIter' -- maximum number of iterations

'StopTolerance' -- stopping tolerance

*Outputs*

x -- last iterate (hopefully an approximate solution)

Out.iter -- # of iterations taken

### GOMP.m

Group Orthogonal Matching Pursuit for l_20 minimization problem:

minimiza ||x||_{2,1}
subject to y = Ax,

*Required inputs*

A -- an m x n matrix

y -- an m x 1 vector

group -- an n-entry vector whose i-th entry is the group number of x_i,
         the group index should start with 1, end with c, where c is the 
         # of group

*Optional inputs*

'maxIter' -- maximum number of iterations

'StopTolerance' -- stopping tolerance

*Outputs*

x -- last iterate (hopefully an approximate solution)

Out.iter -- # of iterations taken

### BOMP.m

Block Orthogonal Matching Pursuit for l_20 minimization problem:

minimiza ||x||_{2,1}
subject to y = Ax,

*Required inputs*

A -- an m x n matrix

y -- an m x 1 vector

group -- an n-entry vector whose i-th entry is the group number of x_i

*Optional inputs*

'maxIter' -- maximum number of iterations

'StopTolerance' -- stopping tolerance

*Outputs*

x -- last iterate (hopefully an approximate solution)

Out.iter -- # of iterations taken
