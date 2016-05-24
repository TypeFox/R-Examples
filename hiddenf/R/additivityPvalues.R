additivityPvalues <-
function(ymtx.out)
{
# this package computes pvalues for tests of additivity using the methods of
###############################################
# Franck and Osborne (2013)                   # 
# Kharrati-Kopaei and Sadooghi-Alvandi (2007) #
# Malik (2014)                                #
###############################################
# Tukey (1949)                                #
# Mandel (1961)                               #
###############################################
# Notes: KKSA does not allow singleton groups #
#        Mandel requires more than two levels #
#        for each factor                      #
###############################################
#ymtx.out <- HiddenF(ymtx)
Mandel.pvalue <- round(MandelPvalue(ymtx.out)$pvalue,digits=4)
Tukey.pvalue <- round(TukeyPvalue(ymtx.out)$pvalue,digits=4)
KKSA.pvalue <- round(KKSAPvalue(ymtx.out)$pvalue,digits=4)
Malik.pvalue <- round(MalikPvalue(ymtx.out)$pvalue,digits=4)
list(
Malik.pvalue = Malik.pvalue, 
Mandel.pvalue =Mandel.pvalue,
Tukey.pvalue =Tukey.pvalue ,
KKSA.pvalue = KKSA.pvalue, 
# hf.pvalue = round(ymtx.out$adjpvalue,digits=4))
ACMIF.pvalue = round(ymtx.out$adjpvalue,digits=4))
}
