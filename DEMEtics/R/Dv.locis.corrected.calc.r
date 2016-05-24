# Function called within p.value.correcture.r:

Dv.locis.corrected.calc <- function(Dv.locis.correction.onelocus){
# For every locus the correcture is carried out separately .
p=as.numeric(as.vector(Dv.locis.correction.onelocus$p.values))
# A vector with the p-values obtained for one locus is extracted.

p.bonferroni <- p.adjust(p,method="bonferroni")
# Bonferroni correction ("bonferroni") in which the p-values are
# multiplied by the number of comparisons.
          
p.holm <- p.adjust(p,method="holm")
# There is no reason to use the unmodified Bonferroni correction
# because it is dominated by Holm's method, which is also valid under
# arbitrary assumptions. It is designed to give strong control of the family
# wise error rate. It is less conservative than the bonferroni correction.
                          
p.hommel <- p.adjust(p,method="hommel")
# Hommel's method is valid when the hypothesis tests are independent
# or when they are non-negatively associated (Sarkar, 1998; Sarkar and
# Chang, 1997). Designed to give strong control of the family wise error rate.
                          
pBH <- p.adjust(p,method="BH")
# method of Benjamini and Hochberg.

Dv.locis.corrected <- cbind(Dv.locis.correction.onelocus,p.bonferroni,p.holm,p.hommel,pBH)
invisible(Dv.locis.corrected)
  
}

