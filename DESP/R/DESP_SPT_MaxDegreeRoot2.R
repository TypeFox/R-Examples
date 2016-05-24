# DESP/R/DESP_SPT_MaxDegreeRoot2.R by A. S. Dalalyan and S. Balmand  Copyright (C) 2015-
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

DESP_SPT_MaxDegreeRoot2 <-
function(Graph) {
  # shortest-path trees - of maximum height equal to 1 - computation choosing root a priori as the node of maximal degree
  
  sptList = list();

  nodesG = nodes(Graph);
  deg = degree(Graph);

  unUsedNodes = c(1:length(nodesG));

  while(length(unUsedNodes) != 0){

    indRoot = unUsedNodes[which.max(deg[unUsedNodes])];
    root = nodesG[indRoot];

    spt = dijkstra.sp(Graph,start=root, eW=unlist(edgeWeights(Graph)));

    unUsedNodes = intersect(unUsedNodes,which(spt$distances == Inf));

    b=(spt$distances != Inf);

    x=unname(spt$penult);

    nsetTree = (x==indRoot | x[x]==indRoot);
    ccTree = as.numeric(names(spt$penult[b]));
    spTree = list(root=indRoot,penult=x,cc=ccTree,nodes=which(nsetTree));
    sptList = append(sptList,list(spTree));

    cut=which(!nsetTree & spt$distances != Inf);
    while (length(cut) > 0){
      grCut = subGraph(snodes=nodesG[cut],graph=Graph);
      newRoot = names(which.max(degree(grCut)));
      sptNew = dijkstra.sp(Graph,start=newRoot, eW=unlist(edgeWeights(Graph)));
      xNew=unname(sptNew$penult);

      nsetTree = (xNew==as.numeric(newRoot) | xNew[xNew]==as.numeric(newRoot));
      nsetTree[-cut] = FALSE;
      spTree = list(root=as.numeric(newRoot),penult=xNew,cc=ccTree,nodes=which(nsetTree));
      sptList = append(sptList,list(spTree));

      cutNew = which(!nsetTree);
      cut = intersect(cut,cutNew);
    }

  }

  return(sptList)
}
