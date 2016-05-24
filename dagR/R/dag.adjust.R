dag.adjust <-
function(dag,A=c())
{ # wrapper for dag.adjustment() and/or find.paths() and eval.paths();
  if(length(A)>0)
  { dag<-dag.adjustment(dag,A);
    dag$searchType <- NULL;
    dag$searchRes <- NULL;
  } else
  { writeLines('The adjustment set is empty. Function dag.adjust does not call function dag.adjustment, but only find.paths and eval.paths.');
  }
  dag<-find.paths(dag);
  dag<-eval.paths(dag);
  return(dag);
}

