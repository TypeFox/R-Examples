library(ergm.userterms)
opttest({
rm(list=ls())
{
data(florentine)
summary(flomarriage~mindegree(3))
summary(flomarriage~mindegree(1,by="priorates"))
fit <- ergm(flomarriage~edges+mindegree(1,by="priorates"))
summary(fit)
}
},"mindegree-term.Rd")
