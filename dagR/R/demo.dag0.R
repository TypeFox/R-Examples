demo.dag0 <-
function()
{ # creates a DAG that was used during development;
  dag<-dag.init(covs=c(1,1), arcs=c(0,2, 1,2, 1,0, -1,2));
  return(dag);
}

