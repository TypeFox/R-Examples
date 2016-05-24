sim.Dvine <-function(DV,N=NULL) {
  
  Dv <- list()
  SSi <- list()
  S <- 1:length(DV$S)
  if(is.null(N)) N <- DV$N
  x <- matrix(runif(N*length(S)),ncol=length(S))
  base <- get("base",DV)
  K <- DV$K
  
  for ( i in 1:length(S))
    {
      vine.knot <-  list(j1=S[i],j2=NULL,D=NULL,v=NULL, u = x[,i] )
      SSi <- append(SSi, list(vine.knot))
    }
  Dv <- append(Dv,list(SSi))
  
  level <- 1
  
  SS <-  Dv[[level]]
  nSS <- length(SS)
  while (nSS>1)
    {
      #print(level)
      SSi <-  list()
      for (i in 2:nSS)
        {
          S1 <- c(SS[[i-1]]$j1, SS[[i-1]]$j2,SS[[i-1]]$D)
          S2 <- c(SS[[i]]$j1, SS[[i]]$j2,SS[[i]]$D)
          index1 <- rep(TRUE,length(S1))
          index2 <- rep(TRUE,length(S2))
          S3 <- c()
          for (j in 1:length(S1)) {
              indexi <- S1[j]==S2
              if (sum(indexi)>0)
                {
                  S3 <- c(S3, S1[j])
                  index1[j] <- FALSE
                  index2[indexi] <- FALSE                
                }
            }
          if(!is.null(S3)) S3 <- sort(S3)
          S1 <- S1[index1]
          S2 <- S2[index2]
                                        #vine.knot <-  list(j1=S1,j2=S2,D=S3,v=NULL,u =NULL)
          vine.knot <- list(j1=S1,j2=S2,D=S3,v=DV[[1]][[level+1]][[i-1]]$v,u=NULL)
                                        #SSi <- list(SSi, list(knoti))
          SSi <- append(SSi, list(vine.knot))
        }
      Dv <- append(Dv, list(SSi))
      level <- level + 1
      SS <-  Dv[[level]]
      nSS <- length(SS)
    }
  Index.basis.D <- matrix(NA,(K+1)^2,2)
  Index.basis.D[,1] <- rep(seq(1,(K+1)),(K+1))
  Index.basis.D[,2] <- sort(Index.basis.D[,1])
 
  x.seq <- seq(0,1,length=N)
  index.b <- matrix(seq(1,N))
  for(j in 2:length(S)) {
    #print(j)
    j.seq <- seq(1,j)
    l.seq <- seq(j,1)
    
    for(l in 2:j) {
      Dv[[j.seq[1]]][[l.seq[1]]]$u <- apply(index.b,1,help.func,j=j,l=l,Dv=Dv,x.seq=x.seq,l.seq=l.seq,j.seq=j.seq,K=K,N=N,Index.basis.D,base=base,q=2)
    }
    for(l in 2:j) {
      Dv[[j.seq[l]]][[l.seq[l]]]$u <- as.vector(cond.cop(cbind(Dv[[j.seq[l-1]]][[l.seq[l]]]$u,Dv[[j.seq[l-1]]][[l.seq[l-1]]]$u),Dv[[j.seq[l]]][[l.seq[l]]]$v,K=K,diff="u2",Index.basis.D,base=base,q=2))
    }
  }
  
  fit <- Dv[[1]][[1]]$u
  for(i in 2:length(S)) fit <- cbind(fit,Dv[[1]][[i]]$u)
  return(fit)
}
