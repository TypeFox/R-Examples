assoc.exists <-
function(dag, a, b)
{ # internally used by dag.adjustment;
  # checks, if an association already exists, i.e. doesn't need to be introduced;
  rv<-FALSE;
  for(i in 1:length(dag$arc[,2]))
  {
    if( ( (dag$arc[i,1]==a && dag$arc[i,2]==b) ||
          (dag$arc[i,2]==a && dag$arc[i,1]==b) 
        ) &&
        (dag$arc.type[i]==1) )
    {
      rv<-TRUE;
    }
  }
  return(rv);
}

