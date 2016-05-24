demo.dag6 <-
function()
{  # creates a miscellaneous DAG;
  dag<-dag.init(covs=c(2,1,1,1,1), arcs=c(1,0, 1,2, 3,2, 3,-1, 4,3, 5,3));
  dag$x<-c(0.000, 0.211, 0.492, 0.492, 0.236, 0.098, 1.000);
  dag$y<-c(0.000, 0.300, 0.300, 0.663, 0.550, 0.816, 0.000);
  return(dag);
}

