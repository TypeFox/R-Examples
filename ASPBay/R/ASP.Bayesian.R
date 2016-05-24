ASP.Bayesian <-
function(N, Tem_Gen, Index_Gen, IBD, snp, thin = 1, sd.freq=0.05, sd.psi=0.05, p0=c(rep(1/4,4),1), psi.prior=0)
{
  Index_Gen <- data.frame(Index_Gen)
  Tem_Gen <- data.frame(Tem_Gen)
  
  if (is.character(snp)) snp <- which(colnames(Tem_Gen)==snp)
  
  S <- c( sum(Tem_Gen[,snp]==0, na.rm=TRUE), sum(Tem_Gen[,snp]==1, na.rm=TRUE), sum(Tem_Gen[,snp]==2, na.rm=TRUE) )
  R <- c( sum(Index_Gen[,snp]==0 & IBD==0, na.rm=TRUE), sum(Index_Gen[,snp]==0 & IBD==1, na.rm=TRUE), sum(Index_Gen[,snp]==0 & IBD==2, na.rm=TRUE),
          sum(Index_Gen[,snp]==1 & IBD==0, na.rm=TRUE), sum(Index_Gen[,snp]==1 & IBD==1, na.rm=TRUE), sum(Index_Gen[,snp]==1 & IBD==2, na.rm=TRUE),
          sum(Index_Gen[,snp]==2 & IBD==0, na.rm=TRUE), sum(Index_Gen[,snp]==2 & IBD==1, na.rm=TRUE), sum(Index_Gen[,snp]==2 & IBD==2, na.rm=TRUE) )

  X <- .Call('ASPBay_MHcpp', PACKAGE = 'ASPBay', N, thin, S, R, sd.freq, sd.psi, p0, psi.prior)
  
  return(list(f_ab=X[,1], f_Ab=X[,2], f_aB=X[,3], f_AB=X[,4], OR=X[,5]))
}
