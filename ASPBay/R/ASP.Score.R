ASP.Score <-
function(Tem_Gen, Index_Gen, IBD)
{
  Index_Gen <- data.frame(Index_Gen)
  Tem_Gen <- data.frame(Tem_Gen)
  
  if (ncol(Tem_Gen)!=ncol(Index_Gen)) return('Controls and cases have not the same number of SNPs')
  
  p <- ncol(Tem_Gen)
  score <- rep(0,p)
  pvalue <- rep(0,p)
  for (i in 1:p)
  {
    S <- c( sum(Tem_Gen[,i]==0, na.rm=TRUE), sum(Tem_Gen[,i]==1, na.rm=TRUE), sum(Tem_Gen[,i]==2, na.rm=TRUE) )
	R <- c( sum(Index_Gen[,i]==0 & IBD==0, na.rm=TRUE), sum(Index_Gen[,i]==0 & IBD==1, na.rm=TRUE), sum(Index_Gen[,i]==0 & IBD==2, na.rm=TRUE),
            sum(Index_Gen[,i]==1 & IBD==0, na.rm=TRUE), sum(Index_Gen[,i]==1 & IBD==1, na.rm=TRUE), sum(Index_Gen[,i]==1 & IBD==2, na.rm=TRUE),
            sum(Index_Gen[,i]==2 & IBD==0, na.rm=TRUE), sum(Index_Gen[,i]==2 & IBD==1, na.rm=TRUE), sum(Index_Gen[,i]==2 & IBD==2, na.rm=TRUE) )

    q_hat <- ( sum(R[4:6]) + 2*sum(R[7:9]) + S[2] + 2*S[3] )/( 2*(sum(R)+sum(S)) )
    U1 <- -4*R[3] - 3*R[2] - 2*R[1] +
          -4*R[6] - 3*R[5] - 2*R[4] +
          -4*R[9] - 3*R[8] - 2*R[7]
    U0 <- 2*R[6] + 3/2*R[5] + R[4] +
          4*R[9] + 3*R[8]   + 2*R[7]
    U <- U1*q_hat + U0
  
    sigma2 <- (1-q_hat)*q_hat*sum(R)*(19*sum(S)+sum(R)-1)/(4*(sum(R)+sum(S)))
	
	score[i] <- U/sqrt(sigma2)
	pvalue[i] <- 2*(1-pnorm(abs(U/sqrt(sigma2))))
  } 
  
  names(score) <- colnames(Tem_Gen)
  names(pvalue) <- colnames(Tem_Gen)
  
  return( list(Value=score, Pvalue=pvalue ) )
}

