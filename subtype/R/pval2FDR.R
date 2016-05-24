pval2FDR <-
function(pval, lim = 0.7)
{
 n= length(pval)
 Fp = rank(pval)/length(pval)
 p0 = sum(pval>lim)/((1-lim)*n)
 p0 = min(p0, 1)

 FDRp = p0 * pmin(pval/Fp, 1)
 ord = order(pval)
 FDR.o = FDRp[ord]
 b = rev(cummin(rev(FDR.o)))
 FDR = rep(0, n)
 FDR[ord] = b
 return(FDR)
}
