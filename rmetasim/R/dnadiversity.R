
                                        #
#These files implement some mismatch distribution calculations.  Incomplete of course.
#base differences  matrix of the sequences in the landscape
#
basediff<-function(lnum=1,Rland)
  {
    if (is.landscape(Rland))
      if (Rland$intparam$locusnum>=lnum)
        if (Rland$loci[[lnum]]$type==253)
          {
            sl<-landscape.locus.states(Rland,lnum);
            rmat<-matrix(0,nrow=length(sl[[1]]),ncol=length(sl[[1]]));
            for (i in 1:length(sl[[1]]))
              for (j in i:length(sl[[1]]))
                {
                  if (i!=j)
                    {
                      vi<-strsplit(sl$state[[i]],NULL)[[1]]
                      vj<-strsplit(sl$state[[j]],NULL)[[1]]
                      rmat[j,i]<-length(vi)-sum(vi==vj);
                      rmat[i,j]<-rmat[j,i];
                    }
                }
            list(sl[[1]],rmat);
          }
  }
#
# produce a table of mismatches for a particular locus
#
landscape.mismatchdist<-function(Rland,lnum=1)
  {
    bd<-basediff(lnum,Rland);
    sl<-bd[[1]];
    dmat<-bd[[2]];
    lt<-landscape.locus(Rland,lnum);
    itbl<-table(lt[,(landscape.democol()+1):ncol(lt)]);
    ttbl<-as.table(table(c(0,seq(max(dmat))))*0);
    for (n in names(itbl))
      {
#        print(paste("Working on: ",n))
        mtbl<-as.table(table(dmat[seq(along=sl)[sl==as.numeric(n)],])*itbl[[n]]);
        for (cn in names(mtbl))
          {
            ttbl[[cn]]<-ttbl[[cn]]+mtbl[[cn]];
          }
      }
    ttbl
  }
