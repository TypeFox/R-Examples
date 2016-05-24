"Q2vect" <-
function(XtX,xbar,E.S,E0l,E1,E2,MT0,MS0,l,u,q1,q2,beta,q4,q5) {
i1    <- integrate(s1phi.w, lower=l,upper=u)$value
i2    <- integrate(s1s2phi.w, lower=l,upper=u)$value
Ex    <- as.matrix(xbar)
Q2    <- Ex%*%(i2+q4*t(E1)%*%t(MT0)%*%Ex+q5*MS0%*%E1-beta*i1)
tmp   <- E2 + q4*E.S%*%t(MT0)%*%Ex +q5*E.S%*%t(MS0) - beta*E0l
Q2    <- Q2 + q1*XtX%*%MT0%*%tmp + q2*Ex%*%MS0%*%tmp
Q2}

