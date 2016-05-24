plotcompprm <-
function(prmdcvobj,...)
{
# Generate plot showing optimal number of components for
# Repeated Double Cross Validation of PRM
#

prm.dcv.optcomp <- table(prmdcvobj$optco)/sum(table(prmdcvobj$optco))
plot(names(prm.dcv.optcomp),prm.dcv.optcomp,type="b",
  xlab="Number of components",ylab="Relative frequency for optimal number",
  cex.lab=1.2,...)
optcomp <- as.numeric(names(which.max(prm.dcv.optcomp)))
abline(v=optcomp,lty=2)

list(optcomp=optcomp,compdistrib=prm.dcv.optcomp)
}

