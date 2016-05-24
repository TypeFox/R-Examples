"Q1mat" <-
function(XtX,xbar,E.S,E1,MT0,MS0,l,u,q1,q2) {
i1    <- integrate(s1s1phi.w, lower=l,upper=u)$value
Ex    <- as.matrix(xbar)
Q1    <- i1*XtX + q1*Ex%*%t(E1)%*%t(MT0)%*%XtX + q2*Ex%*%t(E1)%*%t(MS0)%*%t(Ex)
tmp   <- E1%*%t(Ex) + q1*E.S%*%t(MT0)%*%XtX + q2*E.S%*%t(MS0)%*%t(Ex)
Q1    <- Q1 + q1*XtX%*%MT0%*%tmp
tmp   <- E1%*%t(Ex) + q1*E.S%*%t(MT0)%*%XtX + q2*E.S%*%t(MS0)%*%t(Ex)
Q1    <- Q1 + q2*Ex%*%MS0%*%tmp
Q1}

