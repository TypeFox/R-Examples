demo.dag5 <-
function()
{ # creates a miscellaneous DAG;
  # check out adjustment for the outcome's child!
  dag<-dag.init(covs=c(2,1), arcs=c(1,-1, -1,2));
  return(dag);
}

