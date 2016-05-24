MSLT.S <-
function (rates)  # replaced by M1 destin = 3th dimension
{ numstates <- dim(rates)[2]
  P<- apply(rates,1,function(x) 
	  MatrixExp(x,t=1,n=5,k=3,method="series"))
  zz <- array (P,c(numstates,numstates,nrow(rates)))
  P2 <- aperm(zz,c(3,1,2))  # P2 Age,dest,origin
 dimnames(P2) <- dimnames(rates)
  # S has the transition probabilities from each origin state to the different
  # destination states in conseciutive rows.
  # Ages are in columns
  # The following code rearranges the probabilities 
 S <- array (0,dim=c(nrow(rates),numstates,numstates))
 dimnames (S) = dimnames(rates)  # S Age,dest,orig
 S[1,,] <- diag(1,nrow=numstates,ncol=numstates)
 P3 <- P2   # aperm (P2,c(1,3,2))
 for (ix in 1:(nrow(rates)-1))
  { S[ix+1,,] <- P3[ix,,] %*% S[ix,,]
  }
 # apply(S,c(1,3),sum)  sum = 1 
  class(S) <- "MSLT.S"
  return(list (S=S,
               P=P3))
}
