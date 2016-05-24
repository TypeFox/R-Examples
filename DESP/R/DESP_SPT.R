# DESP/R/DESP_SPT.R by A. S. Dalalyan and S. Balmand  Copyright (C) 2015-
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

DESP_SPT <-
function(X,B,method='1',Theta=NULL) {
  # estimation of the diagonal of the precision matrix using shortest path trees, when the true value of B is known or has already been estimated
  # the observations of the data matrix X are assumed to have zero mean
  # main function

  # read the sample size and the number of variables
  D = dim(X);
  n = D[1];               # n is the sample size
  p = D[2];               # p is the dimension

  # compute the sample cov matrix
  if(is.null(Theta))
    {
    S = crossprod(X)/n;
    }
  else
    {
    S = crossprod(X - Theta %*% MASS::ginv(B))/n;
    }
  
  if (method=='1'){ # considers only the presence or absence of edges to build the shortest path trees
    G = DESP_Weighted_Graph(matrix(1,p,p)*sign(B),n);
    trees = DESP_SPT_MaxDegreeRoot(G);
    Phi = DESP_SPT_Phi(S,B,trees);
  }
  else if (method=='2'){ # chooses the root as the node of maximal degree 
    G = DESP_Weighted_Graph(B,n);
    trees = DESP_SPT_MaxDegreeRoot(G);
    Phi = DESP_SPT_Phi(S,B,trees);
  }  
  else if (method=='2.1'){ # chooses the root as the node of maximal degree and limits the height of shortest path trees to 1
    G = DESP_Weighted_Graph(B,n);
    trees = DESP_SPT_MaxDegreeRoot2(G);
    Phi = DESP_SPT_Phi2(S,B,trees);
  }
  else if (method=='3'){ # get the maximum weighted tree among all shortest path trees for each connected component
    G = DESP_Weighted_Graph(B,n);
    trees = DESP_SPT_MaxWeight(G);
    Phi = DESP_SPT_Phi(S,B,trees);
  }

  return(1/Phi);
}


DESP_SPT_Phi <-
function(S,B,trees) {
  # estimation of the inverse of the diagonal of the precision matrix using shortest path trees

  # read the sample size and the number of variables
  D = dim(S);
  p = D[2];               # p is the dimension

  Phi = c(1:p)*0;

  Dlt = c(1:p)*0;

  penult = trees$penult;
  for (i in c(1:p)){
    j = penult[i];
    m = 1;
    h = i;
    while (j != h){ # j different of the root of the considered tree
      m = m * B[h,j] / B[j,h];
      h = j;
      j = penult[j];
    }
    Dlt[i] = m;
  }

  cc = trees$cc; 
  for(r in c(1:length(cc))){
    dc = Dlt; dc[-cc[[r]]] = 0;
    dr = sum(diag(crossprod(diag(dc),crossprod(S,B))))/(length(cc[[r]]));
    idc = 1/dc; idc[-cc[[r]]] = 0;
    Phi = Phi + dr * idc;
  }

  return(Phi);
}


DESP_SPT_Phi2 <-
function(S,B,trees) {
  # estimation of the inverse of the diagonal of the precision matrix using shortest path trees - of maximum height equal to 1

  # read the sample size and the number of variables
  D = dim(S);
  p = D[2];               # p is the dimension

  Phi = c(1:p)*0;

  for (spt in trees){

    Dlt = c(1:p)*0;

    penult = spt$penult;
    for (i in c(1:p)){
      j = penult[i];
      m = 1;
      h = i;
      while (j != h){ # j different of the root of the considered tree
        m = m * B[h,j] / B[j,h];
        h = j;
        j = penult[j];
      }
      Dlt[i] = m;
    }

    nodes = spt$nodes;
    cc = spt$cc; 
    dc = Dlt; dc[-cc] = 0;
    dr = sum(diag(crossprod(diag(dc),crossprod(S,B))))/(length(cc));
    idc = 1/dc; idc[-nodes] = 0;
    Phi = Phi + dr * idc;

  }

  return(Phi);
}

