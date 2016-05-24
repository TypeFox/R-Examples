"E2vect" <-
function(xbar,kl,ku,l,u) {
i1  <- integrate(s2psiphi.w, lower=kl,upper=ku)$value
i2  <- integrate(s2chiphi.w, lower=l,upper=u)$value
E2  <- matrix(c(i1*xbar,i2),ncol=1)
E2}

