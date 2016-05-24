"gminit" <-
function(aa=c(1,0,1),cc=c(0,0),sigma=NA,alpha=NA,obs) {
if (missing(obs)) {
 if (is.na(sigma)) messagena("obs or sigma")
 if (is.na(alpha)) messagena("obs or alpha")}
else {
 MeanObs <- mean(obs)
 sigma0 <- var(obs)/MeanObs
 alpha0 <- MeanObs/sigma0
 if (is.na(sigma)) sigma <- sigma0
 if (is.na(alpha)) alpha <- alpha0
}
f.res <- .Fortran("gminit",
aa=to.single(aa),
cc=to.single(cc),
sigma=to.single(sigma),
alpha=to.single(alpha))
list(alpha=alpha, sigma=sigma)
}

