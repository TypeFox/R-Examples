# DESP/R/DESP_Weighted_Graph.R by A. S. Dalalyan and S. Balmand  Copyright (C) 2015-
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License (version 3) as published by
#  the Free Software Foundation.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
#

DESP_Weighted_Graph <-
function(B,n) {
  # get a graph representation from the matrix B
  # the resulting graph is a weighted undirected without loops
  # the weights of the edges depend on squared partial correlations

  # read matrix size
  D = dim(B);
  p = D[2];               # p is the dimension

  # compute off-diagonal squared partial correlations
  SPC = DESP_SqPartCorr(B-diag(p),n);

  # matrix of weights computation
  W = exp(-SPC)*(SPC!=0);

  # conversion to a sparse matrice stored in compressed sparse row format
  W = as.matrix.csr(W);

  # graph representation build from a sparse matrice
  # graphNEL : This is a class of graphs that are represented in terms of nodes and an edge list. This is a suitable representation for a graph with a large number of nodes and relatively few edges. 
  G = sparseM2Graph(W, nodeNames=as.character(c(1:p)), edgemode="undirected")

  return(G)
}
