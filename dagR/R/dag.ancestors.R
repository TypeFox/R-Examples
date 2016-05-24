dag.ancestors <-
function(dag,A)
{ # identify ancestors of A;
  ancest<-c();
  rv<-c(A);

  # add the X-Y arrow for determination of ancestors;
  # also keep all other arrows emenating from X;
  arcs<-rbind(dag$arc, c(1,length(dag$names)));
  types<-c(dag$arc.type, 0);
  
  # search for immediate ancestors (parents)
  for(i0 in 1:length(A))
  {
    i<-0; 
    while(i<length(arcs[,2]))
    { i<-i+1;
      # only use directed arrows to identify ancestors;
      if(arcs[i,2]==A[i0] && types[i]!=1)
      { ancest<-c(ancest, arcs[i,1]);
      }
    }
  }
  if(length(ancest)>0)
  { rv<-c(A,dag.ancestors(dag,ancest));
  }
  return(unique(rv));
}

