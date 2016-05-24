cronbach <-
function(v1)
{
v1 <- na.omit(v1)
nv1 <- ncol(v1)
pv1 <- nrow(v1)
alpha <- (nv1/(nv1-1))*(1 - sum(apply(v1,2,var))/var(apply(v1,1,sum)))
resu <- list("sample.size"=pv1,"number.of.items"=nv1,"alpha"=alpha)
resu
}
