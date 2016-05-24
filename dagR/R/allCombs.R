allCombs <-
function(x, force=c(), trace=FALSE)
{ # internally used by brute.search;
  # creates eligible adjustment set combinations;
  combs<-NULL;

  # Remove the variables, that are forced to stay in,
  # from the set from which permutations are created.
  # The forced ones will be put in front of each premutation.
  if(length(force)>0)
  { remove.x<-c();
    for(i in 1:length(x))
    { if(trace==TRUE) writeLines(paste(x[i]));
      if(is.in(as.numeric(x[i]), c=force)==TRUE)
      { remove.x<-c(remove.x, i);
        if(trace==TRUE) writeLines('removing');
    } }
    if(length(remove.x)>0)
    { x<-x[-remove.x];
  } }

  if(is.null(x)==FALSE)
  { for(i in 0:length(x))
    { if(length(x)>1)
      { new.combs<-t(combn(x, i));
      } else
      { if(i==0) { new.combs<-as.matrix(NA);
        } else
        { new.combs<-as.matrix(x);
        }
      }
      new.combs<-cbind(new.combs,
                       matrix(rep(NA, dim(new.combs)[1]*
                                          (length(x)-dim(new.combs)[2])),
                              nrow=dim(new.combs)[1]));
      combs<-rbind(combs, new.combs);
  } }

  if(length(force)!=0)
  {  combs<-cbind(matrix(rep(c(force), dim(combs)[1]),
                         nrow=dim(combs)[1], byrow=TRUE),
                  combs);
  }
  return(combs);
}

