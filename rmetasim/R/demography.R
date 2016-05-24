#
#
#
landscape.demography <- function(Rland)
{
  #this routine is not ready for use.
  rl <- vector("list",length=Rland$intparam$habitats+1)
  A.full <- Rland$demography$epoch[[Rland$intparam$currentepoch+1]]$R + Rland$demography$epoch[[Rland$intparam$currentepoch+1]]$S
  oblanddist <- table(Rland$individuals[,1])
  rlcnt <- 1
  for (i in 0:(Rland$intparam$habitats-1))
    {

      loclist <- NULL
      strt <- (i*Rland$intparam$stages)+1
      len <- Rland$intparam$stages-1
      obdist <- c(oblanddist[(strt):(strt+len)])
#print(paste("i",i,"strt",strt,"len",len))
#      print(obdist)
      A <- NULL
      A <- A.full[(strt):(strt+len),(strt):(strt+len)]
      es <- eigen(A)
      l <- es$values
#      print(A)
      statedist <- es$vectors[,(which(Re(l)==max(Re(l))))]
      statedist <- (statedist/sum(statedist))
#      print(statedist)
      exdist <- statedist*sum(obdist)
      chisq <- sum(((obdist-exdist)^2)/exdist)
      df <- length(obdist)-1
      prob <- 1-pchisq(chisq,df)
      loclist <- list(lambda=l[which(Re(l)==max(Re(l)))],statedist=statedist,
                      obdist=obdist,exdist=exdist,chisq=chisq,df=df,prob=prob)
      rl[[rlcnt]] <- loclist
      rlcnt <- rlcnt+1
    }
  A <- NULL
  es <- NULL
  A <- A.full
  es <- eigen(A)
  l <- es$values
  print(A)
  statedist.1 <- es$vectors[,(which(Re(l)==max(Re(l))))]
#  print(statedist.1)
  statedist.2 <- (statedist.1/sum(statedist.1))
  exdist <- statedist.2*sum(oblanddist)
print("exdist")
  
  print(exdist)
print("oblanddist")
  print(oblanddist)
  chisq <- sum(((oblanddist-exdist)^2)/exdist)
  df <- length(oblanddist)-1
  prob <- 1-pchisq(chisq,df)

  loclist <- list(lambda=l[which(Re(l)==max(Re(l)))],statedist=statedist.2,obdist=oblanddist,exdist=exdist,chisq=chisq,df=df,prob=prob)
  rl[[length(rl)]] <- loclist 
  rl
}
