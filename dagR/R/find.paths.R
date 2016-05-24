find.paths <-
function(dag)
{ # identifies backdoor paths, open and closed ones;
  finish<-FALSE;
  backdoorpaths<-NULL;
  pathsN<-0;
  zaehler<-0;
  cpos<-0;                 # cpos = current position
  X<-1;                    # exposure node
  Y<-length(dag$names);    # outcome node
  
    # search initial backdoor path
    used<-c(TRUE, rep(FALSE, length(dag$names)-1));
    cpos<-X;
    cpath<-c();                           # cpath = current path
    cycledPath<-c(0);  # Stores, how far I cycled through to move on from node.
    finish2<-FALSE;
    if(length(dag$arc[,1])==0)
    { finish2<-TRUE;
    }
    stepback<-FALSE;
    search1<-1;
    while (finish2==FALSE)
    {
      if(     # path not connected.
              (dag$arc[search1,1]!=cpos && dag$arc[search1,2]!=cpos)
              # path connected but target node already used.
          ||  (dag$arc[search1,1]==cpos && used[dag$arc[search1,2]]==TRUE)
          ||  (dag$arc[search1,2]==cpos && used[dag$arc[search1,1]]==TRUE)
              # path arrived at outcome.
          ||  (cpos==Y)
              # already cycled too far.
          ||  (stepback==TRUE)
        )
      {
        # If i'm at X && starting only && already cycled through all arcs, quit.
        if(cpos==X && length(cpath)==0 && search1==length(dag$arc[,1]))
        { finish2<-TRUE;
        } else
        {        # If i have arrived at Y and need to search next path.
                 # OR
                 # If i'm anywhere (e.g. returned to X via a loop),
                 # and can't move, and checked all paths,
                 # go one step back (rechange cpos, rechange used,
                 #                   shorten cpath, shorten cycledPath) and
                 # search from next arc on, only!
          if( cpos==Y || stepback==TRUE || search1==length(dag$arc[,1]) )
          {
            stepback<-FALSE;
            used[cpos]<-FALSE;
            lastStep<-cpath[length(cpath)];
            if(cpos!=dag$arc[lastStep,1])
            {        cpos<-dag$arc[lastStep,1]; 
            } else { cpos<-dag$arc[lastStep,2];
            }
            cpath<-cpath[-length(cpath)];
            search1<-cycledPath[length(cycledPath)-1];
            cycledPath<-cycledPath[-length(cycledPath)];
            cycledPath[length(cycledPath)]<-0;
          }
        }
      } else
      {      # Arc is valid.
             # Add it to the current path.
             # Update the current position.
             # Mark new node as used.
             # Reset the search counter, as all arcs are available again.
             # Store it in cycledPath, and elongate cycledPath.
        cpath<-c(cpath, search1);
        if(cpos!=dag$arc[search1,1])
        {        cpos<-dag$arc[search1,1];
        } else { cpos<-dag$arc[search1,2];
        }
        used[cpos]<-TRUE;
        cycledPath[length(cycledPath)]<-search1;
        cycledPath<-c(cycledPath, 0);
        search1<-0;
        if(cpos==Y)  # If it reaches OUTCOME, write it into backdoorpaths!
        { goodpath<-c(cpath, rep(NA, length(dag$arc[,1])-length(cpath)));
          backdoorpaths<-rbind(backdoorpaths, goodpath, deparse.level=0);
          pathsN<-pathsN+1;
        }
      }
      search1<-search1+1;
      if(search1>length(dag$arc[,1]))
      { search1<-search1-1;
        stepback<-TRUE;
      }
    }

  # delete paths, that start with a directed arc from the exposure
  # (i. e. you do not want to consider exposure effects as backdoor)
  if(is.null(backdoorpaths)==FALSE)
  {
    i<-0;
    del.vec<-rep(TRUE, length(backdoorpaths[,1]));
    while(i < length(backdoorpaths[,1]))  
    { 
      i<-i+1;
      if( (     dag$arc[backdoorpaths[i,1],1] == 1) &&
          (dag$arc.type[backdoorpaths[i,1]]   != 1) )
      {
        del.vec[i]<-FALSE;
        pathsN<-pathsN-1;
      } 
    }
    backdoorpaths<-backdoorpaths[del.vec,];
  }

  dag$pathsN<-pathsN;
  if(dag$pathsN==1){ backdoorpaths<-matrix(backdoorpaths, nrow=1); }

  # some of the later function require an NA as the last value of the paths;
  if(is.null(backdoorpaths)==FALSE)
  { backdoorpaths<-cbind(backdoorpaths, rep(NA, nrow(backdoorpaths)));
  }

  dag$paths<-backdoorpaths;
  
  return(dag);
}

