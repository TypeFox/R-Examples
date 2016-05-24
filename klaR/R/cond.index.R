cond.index <- function(formula, data, ...)
{
reg<-lm(formula, data, x=TRUE, ...)
x.scale <- reg$x%*%diag(1/sqrt(colSums(reg$x^2)))
x.scale.svd<- svd(t(x.scale)%*%x.scale)
ci <- sqrt(max(x.scale.svd$d)/x.scale.svd$d)
return(ci)
}
