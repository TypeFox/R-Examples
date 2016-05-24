demo.dag3 <-
function()
{ # re-creates DAG example 3 in Knueppel,
  # DAG v0.11 documentation Oct 21, 2009;
  dag<-dag.init(covs=c(1,1,1,1), cov.names=c("A","B","C","Z"),
                arcs=c(1,0, 1,4, 2,4, 2,-1, 3,4, 3,-1, 4,0, 4,-1));
  dag$x<-c(0.000, 0.000, 0.501, 1.058, 0.327, 1.000);
  dag$y<-c(0.000, 0.495, 0.574, 0.491, 0.260, 0.000);
  return(dag);
}

