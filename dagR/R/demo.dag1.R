demo.dag1 <-
function()
{ # re-creates DAG in figure 3 of Fleischer & Diez Roux,
  # J Epidemiol Community Health 2008;62:842;
  dag<-dag.init(y.name="incident CVD", x.name="neighborhood violence",
       covs=c(1,1,1),
       cov.names=c("urban residence", "participation", "family history"),
       arcs=c(1,0, 1,2, 3,2, 3,-1));
  dag$x[3]<-0.505; dag$y[3]<-0.269;
  return(dag);
}

