"E1vect" <-
function(xbar,kl,ku,l,u) {
i1  <- integrate(s1psiphi.w, lower=kl,upper=ku)$value
i2  <- integrate(s1chiphi.w, lower=l,upper=u)$value
E1  <- matrix(c(i1*xbar,i2),ncol=1)
E1}

