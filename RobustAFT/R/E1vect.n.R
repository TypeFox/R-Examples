"E1vect.n" <-
function(xbar,ku,u) {
i1  <- integrate(s1psiphi.n, lower=-ku,upper=ku)$value
i2  <- 0 #integrate(s1chiphi.n, lower=-u,upper=u)$value
E1  <- matrix(c(i1*xbar,i2),ncol=1)
E1}

