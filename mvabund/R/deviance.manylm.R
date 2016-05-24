############################################################################
# return the deviance of a manylm object							
############################################################################

deviance.manylm <- function (object, na.action="na.omit", ...) {

if (is.null(object$weighted.residuals))
   wr <- as.matrix(object$residuals)
else
   wr <- as.matrix(object$weighted.residuals)


if (na.action=="na.fail") wr <- na.fail(wr)  else
if (na.action=="na.omit") wr <- na.omit(wr)  else
if (na.action=="na.exclude") wr <- na.exclude(wr) 

# dev <- t(wr) %*% wr 
# return(diag(dev))
dev <- wr^2
drop(rep.int(1, nrow(dev)) %*% dev)

}

