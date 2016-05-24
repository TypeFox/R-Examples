"ntdt.q" <-
function(g,m,alpha=0.00000005,power=0.80)
{
# Creates a table of number of families needed, for
# different values of q, and three different values of
# LD (100%, 75%, 50% of Dmax).
#   ntdt.q(g=1.5,m=0.5,alpha=0.00000005,power=0.80)
 tb <- NULL
 qv <- seq(0.01,0.8,0.01) 
 for (q in qv)
  {
  op <- ntdt(q,m,ld=1,g=g,power=power,alpha=alpha)
  op.75 <- ntdt(q,m,ld=0.75,g=g,power=power,alpha=alpha)
  op.5 <- ntdt(q,m,ld=0.50,g=g,power=power,alpha=alpha)
  tb <- rbind(tb,c(op$nfam,op.75$nfam,op.5$nfam,log10(op$nfam),log10(op.75$nfam),log10(op.5$nfam)))
  }
 tb <- data.frame(cbind(qv,tb))
 colnames(tb) <- c("q","dmax","dmax.75","dmax.50","log.dmax","log.dmax.75","log.dmax.50")
 return(tb)
}

