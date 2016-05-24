"E2vect.n" <-
function(xbar,ku,u) {
i1  <- 0 #integrate(s2psiphi.n, lower=-ku,upper=ku)$value
i2  <- integrate(s2chiphi.n, lower=-u,upper=u)$value
E2  <- matrix(c(i1*xbar,i2),ncol=1)
E2}

