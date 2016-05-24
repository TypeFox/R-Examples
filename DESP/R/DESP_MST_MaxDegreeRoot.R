# DESP/R/DESP_MST_MaxDegreeRoot.R by A. S. Dalalyan and S. Balmand  Copyright (C) 2015-
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

DESP_MST_MaxDegreeRoot <-
function(Graph) {
  # minimum spanning trees computation choosing root as the node of maximal degree
  
  # minimum spanning tree
  r = try(mst <- mstree.kruskal(Graph), silent = TRUE )
  if (inherits (r , "try-error")) {
    print("Warning : error when calling mstree.kruskal in DESP_MST_MaxDegreeRoot : every node is isolated");
    penult = 1:numNodes(Graph);
    cc = penult; for (i in cc){cc[i]=list(i)}
    sptList = list(penult=penult,cc=cc);
    return(sptList)
  }

  dfMst <- data.frame(from=mst$edgeList[1,],to=mst$edgeList[2,],weight=as.vector(mst$weights));
  grMst <- graphBAM(dfMst, edgemode = "undirected");
  treesGraph = as(grMst,"graphNEL");

  # connected components 
  cc = connComp(Graph);
  ccN = list();

  nbNodes = numNodes(Graph);
  penult = c(1:nbNodes);
  usedNodes = vector(mode = "logical", length = nbNodes);
  for (i in 1:length(cc)){
    comp = cc[[i]];
    n = length(comp);
    if (n > 1){
      root = as.numeric(names(which.max(degree(subGraph(snodes=comp,graph=treesGraph)))));
      usedNodes[root] = TRUE;
      fifo = c(root);
      j = 1;
      while (j<n){
        node = fifo[1];
        fifo = fifo[-1];
        for (adj in edges(treesGraph)[[as.character(node)]]){
          adj = as.numeric(adj);
          if (!usedNodes[adj]){
            penult[adj] = node;
            usedNodes[adj] = TRUE;
            fifo = c(fifo,adj);
            j = j + 1;
          }
        }
      }
    }
    ccN = append(ccN,list(as.numeric(comp)));
  }

  sptList = list(penult=penult,cc=ccN);

  return(sptList)
}
