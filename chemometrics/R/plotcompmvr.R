plotcompmvr <-
function(mvrdcvobj,...)
{
# Generate plot showing optimal number of components for
# Repeated Double Cross Validation 
#

pcr.dcv.optcomp=table(mvrdcvobj$optco)/sum(table(mvrdcvobj$optco))
plot(names(pcr.dcv.optcomp),pcr.dcv.optcomp,type="b",
  xlab="Number of components",ylab="Relative frequency for optimal number",
  cex.lab=1.2,...)
optcomp=as.numeric(names(which.max(pcr.dcv.optcomp)))
abline(v=optcomp,lty=2)

list(optcomp=optcomp,compdistrib=pcr.dcv.optcomp)
}

