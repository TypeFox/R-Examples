eval.paths <-
function(dag)
{ # checks if paths are open or blocked;
  # if both blocked "by collider" and "by adjustment", collider dominates;
  if(is.null(dag$paths)==TRUE || dag$pathsN==0)
  { path.status<-NULL;
  } else
  { path.status<-rep('open', dag$pathsN); # first, all are open
    #if(dag$pathsN==1)
    #{ dag$paths<-matrix(dag$paths, nrow=1);
    #}
    for(i in 1:dag$pathsN)
    {
      i2<-1;
      while(is.na(dag$paths[i, i2+1])==FALSE) # continue while another arrow follows
      { if( # collider present if arcs i and i2 point to same node
               (dag$arc[dag$paths[i, i2], 2] == dag$arc[dag$paths[i, i2+1], 2])
            # and neither one is an association
            && (dag$arc.type[dag$paths[i, i2]]   != 1)
            && (dag$arc.type[dag$paths[i, i2+1]] != 1) )
        { path.status[i]<-'blocked by collider';
          i2<-length(dag$paths[i,])-2;
        } else
        { if(   (is.in( as.numeric( dag$arc[dag$paths[i, i2], 1]  ), dag$adj) == TRUE)
             || (is.in( as.numeric( dag$arc[dag$paths[i, i2], 2]  ), dag$adj) == TRUE)
             || (is.in( as.numeric( dag$arc[dag$paths[i, i2+1], 1]), dag$adj) == TRUE)
             || (is.in( as.numeric( dag$arc[dag$paths[i, i2+1], 2]), dag$adj) == TRUE) )
          { 
            path.status[i]<-'blocked by adjustment';
          }
        }
        i2<-i2+1;
      }
    }
  }
  dag$path.status<-path.status;
  return(dag);
}

