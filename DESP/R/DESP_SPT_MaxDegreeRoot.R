# DESP/R/DESP_SPT_MaxDegreeRoot.R by A. S. Dalalyan and S. Balmand  Copyright (C) 2015-
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

DESP_SPT_MaxDegreeRoot <-
function(Graph) {
  # shortest path trees computation choosing root a priori as the node of maximal degree
  
  sptList = list();

  nodesG = nodes(Graph);
  deg = degree(Graph);

  unUsedNodes = c(1:length(nodesG));

  penult = c();
  cc = list();
  while(length(unUsedNodes) != 0){
    indRoot = unUsedNodes[which.max(deg[unUsedNodes])];
    root = nodesG[indRoot];

    spt = dijkstra.sp(Graph,start=root, eW=unlist(edgeWeights(Graph)));

    unUsedNodes = intersect(unUsedNodes,which(spt$distances == Inf));

    b=(spt$distances != Inf);
    spt = spt$penult[b];
    penult = c(penult,spt);
    cc = append(cc,list(as.numeric(names(spt))));
  }
  ord = order(as.numeric(names(penult)));
  sptList = list(penult=unname(penult)[ord],cc=cc);

  return(sptList)
}
