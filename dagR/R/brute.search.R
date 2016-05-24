brute.search <-
function(dag, allow.unknown=FALSE, trace=TRUE, stop=0)
{ # brute search for adjustment sets,
  # i.e. all possible adjustment sets are evaluated;
  # default is that even "unknown" covariables are evaluated;
  # if the input DAG is adjusted, only sets incl. these adjusted covs are evaluated;
  # stop: 0 means no stopping, 1 means stop after 1st sufficient set;
  covsN<-length(dag$names)-2;
  allSets<-c();
  totalPaths<-c();
  openPaths<-c();
  if(covsN>0)
  {
    if(is.null(dag$adj)==TRUE)
    { force<-c();
    } else
    { force<-dag$adj;
    }

    if(allow.unknown==TRUE)
    { allSets<-allCombs(2:(1+covsN), force);
    } else
    { allCovariables<-c(2:(1+covsN));
      theUnknowns<-sapply(X=allCovariables,
                          function(x){ is.unknown(x, dag);})
      theKnowns<-theUnknowns==FALSE;
      knownCovs<-allCovariables[theKnowns];
      allSets<-allCombs(knownCovs, force);
    }

    for(i in 1:dim(allSets)[1])
    { currentSet<-as.vector(na.omit(allSets[i,]));
      dag.temp<-dag.adjust(dag, currentSet);

      totalPaths<-c(totalPaths, dag.temp$pathsN);

      if(dag.temp$pathsN==0)
      { openPaths<-c(openPaths, 0);
      } else
      { openPaths<-c(openPaths,
                     dag.temp$pathsN-sum(dag.temp$path.status!='open'));
      }

      if(trace==TRUE)
      { writeLines(paste(c('set:', currentSet), collapse=' '));
        writeLines(paste(c('paths:', dag.temp$pathsN,
                           '; open:', openPaths[length(openPaths)]),
                         collapse=' '));
      }

      if( openPaths[length(openPaths)]==0 && stop==1 )
      { totalPaths<-c(totalPaths, rep(NA, dim(allSets)[1]-length(openPaths)));
        openPaths<-c(openPaths, rep(NA, dim(allSets)[1]-length(openPaths)));
        break();
      }
    }
  }
  rv<-data.frame(allSets, totalPaths, openPaths);
}

