MSLT.e <-
function (SS,radix)
{ 	S<-SS$S
	nage <- nrow(S)
   iradix <- which(radix>0)
   lt <- apply(S[,,iradix],1,sum) # total survival prob: starts in iradix (1)
    namstates <- unlist(unname(dimnames(S)[2]))
    numstates <- length (namstates)
   LL <- array(0,c(nage,numstates,numstates))
   dimnames(LL) <- dimnames(S)
   for (ix in 1:(nage-1))
   { LL[ix,,] <-  0.5 * (S[ix,,]+S[ix+1,,])
   }# Expected sojourn time beyond age 0
    e0 <- apply(LL,c(2,3),sum)
 #  population based: average person aged 50
 #   e50.p <- apply (LL[51:nage,,],c(2,3),sum)/lt[51]
 # status-based: by status at age 50
 #   e50.s <- apply (LL[51:nage,,],c(2,3),sum)%*%solve(S[51,,])
    e.p <- array(NA,dim=c(nage,numstates,numstates),dimnames=dimnames(S))
    e.s <- e.p
    for (ix in 1:(nage-1))
    {  e.p[ix,,] <- apply (LL[ix:nage,,],c(2,3),sum)/lt[ix]}
    
    LLL <- LL
    for (ix in 1:(nage-1))
    { S[ix,,] <- 0 
      diag(S[ix,,]) <- 1
      for (iy in ix:(nage-1))
      {  S[iy+1,,] <- SS$P[iy,,]%*%S[iy,,]
      	 LLL[iy,,] <- 0.5 * (S[iy,,]+S[iy+1,,])
      }
      e.s[ix,,] <- apply(LLL[ix:nage,,],c(2,3),sum)
    }
     #  zx <- det(S[ix,,])
     #  if (zx < 1e-16) {e.s[ix,,] <- NA} else
     #    {e.s[ix,,] <- apply (LL[ix:nage,,],c(2,3),sum)%*%solve(S[ix,,])}
    
  return (list(L = LL,
               e0=e0, 
               e.p = e.p,
               e.s=e.s ))
}
