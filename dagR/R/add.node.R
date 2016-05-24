add.node <-
function(dag, name="unknown", type=1, x=NA, y=NA)
{ # this adds a node in a convenient way;
  # it requires a function because the node is inserted
  # at the 2nd-to-last position in the vectors...;
  nodesN<-length(dag$x);
  if(is.na(x)) x<-0.5*(min(dag$x)+max(dag$x))+sin(nodesN)*(max(dag$x)-min(dag$x))/4;
  if(is.na(y)) y<-0.5*(min(dag$y)+max(dag$y))+cos(nodesN)*(max(dag$y)-min(dag$y))/4;
  dag$cov.types<-c(dag$cov.types[1:(nodesN-1)], type, dag$cov.types[nodesN]);
  dag$x<-c(dag$x[1:(nodesN-1)], x, dag$x[nodesN]);
  dag$y<-c(dag$y[1:(nodesN-1)], y, dag$y[nodesN]);
  dag$names<-c(dag$names[1:(nodesN-1)], name, dag$names[nodesN]);
  dag$arc<-matrix(sapply(X=dag$arc,
                  FUN=function(x){if(x==nodesN) nodesN+1 else x;}),
                  byrow=FALSE, ncol=2);
  return(dag);
}

