#
# formats a landscape locus to calculate amova
#
landscape.amova.locus <- function(rland,l=1)
{
  loc <- landscape.locus(rland,l)
  if (landscape.ploidy(rland)[l]>1)
    {
      sortloc <- t(apply(loc[,(landscape.democol()+1):(landscape.democol()+2)],1,sort))
      df <- data.frame(matrix(table(paste(sortloc[,1],sortloc[,2]),landscape.populations(rland)),ncol=length(unique(landscape.populations(rland))),byrow=F))
      names(df) <- colnames(table(paste(sortloc[,1],sortloc[,2]),landscape.populations(rland)))
      rownames(df) <- rownames(table(paste(sortloc[,1],sortloc[,2]),landscape.populations(rland)))
    }
  else
  {
    df <- data.frame(matrix(table(paste(loc[,landscape.democol()+1]),landscape.populations(rland)),
                            ncol=length(unique(landscape.populations(rland))),byrow=F))
    names(df) <- colnames(table(paste(loc[,landscape.democol()+1]),landscape.populations(rland)))
    rownames(df) <-   rownames(as.matrix(table(paste(loc[,landscape.democol()+1]),landscape.populations(rland))))
  }
  amova(samples=df)
}


landscape.amova.pairwise <- function (rland) 
{
  habs <- unique(landscape.populations(rland))
  if (length(habs) > 0)
    {
      retmat <- matrix(NA, nrow = length(habs), ncol = length(habs))
      for (i in 1:length(habs))
        for (j in 1:i - 1)
          {
            compland <- rland
            compland$individuals <- compland$individuals[landscape.populations(compland) %in% c(habs[i], habs[j]), ]
            tmpvec <- rep(NA, length(rland$loci))
            n <- 1
            for (l in 1:length(rland$loci))
              {
                tmpvec[l] <- landscape.amova.locus(compland, l)$statphi
                if (!is.na(tmpvec[l])) 
                  n <- n + 1
              }
            retmat[i, j] <- sum(unlist(tmpvec), na.rm = T)/n
          }
      colnames(retmat) <- habs
      rownames(retmat) <- habs
      retmat
    } else {warning("no occupied habitats"); NULL}

}

landscape.amova.pairwise.old <- function(rland)
  {
    retmat <- matrix(NA,nrow=rland$intparam$habitat,ncol=rland$intparam$habitat)
    for (i in 1:rland$intparam$habitat)
      for (j in 1:i-1)
        {
          if (length(landscape.populations(rland))>0)
            {
              if ((max(landscape.populations(rland)==i)>0)&(max(landscape.populations(rland)==j)>0))
                {
                  compland <- rland
                  compland$individuals <- compland$individuals[landscape.populations(compland) %in% c(i,j),]
                  tmpvec <- rep(NA,length(rland$loci))
                  n <- 1
                  for (l in 1:length(rland$loci))
                    {
                      tmpvec[l] <- landscape.amova.locus(compland, l)$statphi
                      if (!is.na(tmpvec[l]))
                        n <- n+1;
                    }
             
                  retmat[i,j] <- sum(unlist(tmpvec),na.rm=T)/n
                }
            }
        }
    retmat
  }


landscape.amova <- function(rland,np=NULL,ns=NULL) #np is the sample size per population.  Only pops with actual > np
  {
    unlist(sapply(1:length(rland$loci),function(x,l){landscape.amova.locus(l,x)$statphi
                                            },l=landscape.sample(rland,np=np,ns=ns)))
  }

