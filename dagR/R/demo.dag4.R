demo.dag4 <-
function()
{ # creates a miscellaneous DAG;
  # check out adjustment for the exposure's child!
  dag<-dag.init(covs=c(1,1,1), arcs=c(0,1, 0,2, 1,2, 1,3, 3,-1));
  dag$y[3]<-0.25;
  return(dag);
}

