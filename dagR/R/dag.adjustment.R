dag.adjustment <-
function(dag,A=NULL)
{ # called internally by dag.adjust();
  # identifies the associations introduced by adjusting for A;

  dag$adj<-A;

  dag$pathsN<-NULL;
  dag$paths<-NULL;

# find all ancestors of A
  ancs<-dag.ancestors(dag, A);

# put in all associations for all ancestors
  # for each ancestor
  for(i in 1:length(ancs))
  {
    # 1. Look up all parents.
    #     For this step, I also have to include the x->y arrow,
    #     as is also done inside dag.ancestors().
    arcs<-rbind(dag$arc, c(1,length(dag$names)));
    types<-c(dag$arc.type, 0);
    parents<-c();
    i2<-0; 
    while(i2<length(arcs[,2]))
    { i2<-i2+1;
      if(arcs[i2,2]==ancs[i] && types[i2]!=1)
      { parents<-c(parents, arcs[i2,1]);
      }
    }
    # 2. if >1 parents, introduce associations if none yet present
    if(length(parents)>1)
    { # examine--sort of--a triangular matrix of the parents
      for(i3 in 2:length(parents))
      { for(i4 in 1:(i3-1))
        { # is there already an association? 
          if(assoc.exists(dag, parents[i3], parents[i4])==FALSE)
          {
            dag$arc<-rbind(dag$arc, c(parents[i3], parents[i4]));
            dag$arc.type<-c(dag$arc.type, 1);
          }
        }
      }
    }    
  }
  return(dag);
}

