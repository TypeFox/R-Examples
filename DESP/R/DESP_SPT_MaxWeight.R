# DESP/R/DESP_SPT_MaxWeight.R by A. S. Dalalyan and S. Balmand  Copyright (C) 2015-
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

DESP_SPT_MaxWeight <-
function(Graph) {
  # maximum weighted tree among all shortest path trees computation 
  
  sptList = list();
  tmpList = list();

  nodesG = nodes(Graph);
  unUsedNodes = c(1:length(nodesG));

  maxLength = unUsedNodes*NA;

  root = NULL;
  minLength = Inf;
  for (i in nodesG){
    spt = dijkstra.sp(Graph,start=i, eW=unlist(edgeWeights(Graph)));
    tmpList = append(tmpList,list(spt));
    mL = max(spt$distances[spt$distances!=Inf]);
    maxLength[as.numeric(i)] = mL;
    if (mL < minLength){
      minLength = mL;
      root = as.numeric(i);
    }
  }

  penult = c();
  cc = list();
  while(length(unUsedNodes) != 0){
    spt = tmpList[[root]];

    unUsedNodes = intersect(unUsedNodes,which(spt$distances == Inf));

    b=(spt$distances != Inf);
    spt = spt$penult[b];
    penult = c(penult,spt);
    cc = append(cc,list(as.numeric(names(spt))));

    root = unUsedNodes[which.min(maxLength[unUsedNodes])];
  }

  ord = order(as.numeric(names(penult)));
  sptList = list(penult=unname(penult)[ord],cc=cc);

  return(sptList)
}
