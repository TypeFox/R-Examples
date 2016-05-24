#=====================================================================================
#
#           Filename:  qqplot.R
#
#        Description:  R function "qqplot.plot"
#
#            Version:  1.0
#            Created:  11-Mar-2009
#           Revision:  none
#  last modification:  15-Jan-2010
#
#             Author:  Maksim V. Struchalin
#        Modified by:  
#            Company:  ErasmusMC, Epidemiology & Biostatistics Department, Department of Forensic Molecular Biology, The Netherlands.
#              Email:  m.struchalin@erasmusmc.nl
#            license:  GPL (>=2)
#
#=====================================================================================

"qqplot.plot" <-
function(chi2, filename="", df=1)
{


#print(cat("Usage: qqplot.plot(chi2, filename)\nInput parametres are simple vector with chi2 filename.\nIf there are no filename then output is on screen.\n\nauthor: M.Struchalin, m.struchalin@erasmusmc.nl\n"))



	chi_observed <- chi2
	chi_observed <- na.omit(chi_observed)
	chi_observed[chi_observed < 0] <- 0	

	chi_expected <- rchisq(length(chi_observed),df=df)
	
	max_chi2_1 <- max(chi_expected, na.rm=T)
	max_chi2_2 <- max(chi_observed, na.rm=T)
	max <- max_chi2_1
	if(max < max_chi2_2)
		{
		max <- max_chi2_2
		}
	max <- round(max+2)

	if(filename != "")
		{
		bitmap(paste(filename, "_qqplot_chi2", ".jpeg", sep=""), type="jpeg")
		}
	qqplot(chi_expected, chi_observed, ylim=c(0,max), xlim=c(0,max), main="qqplot", xlab="chi2 expected", ylab="chi2 observed")
	abline(0, 1,  col="red")
	if(filename != "")
		{
		dev.off()
		}

	if(filename != "")
		{
		bitmap(paste(filename, "_qqplot_pval", ".jpeg", sep=""), type="jpeg")
		}
	pval_expected <- -log10(pchisq(chi_expected, df=df, lower.tail=F))
	pval_observed <- -log10(pchisq(chi_observed, df=df, lower.tail=F))
	pval_max <- -log10(pchisq(max, df=df, lower.tail=F))
	qqplot(pval_expected, pval_observed, ylim=c(0,pval_max), xlim=c(0,pval_max), main="qqplot_pval", xlab="-log10(pval) expected", ylab="-log10(pval) observed")
	abline(0, 1,  col="red")
	if(filename != "")
		{
		dev.off()
		}

}
