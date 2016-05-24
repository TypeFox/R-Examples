"Q3sca.n" <-
function(xbar,E.S,E0l,E2,MT0,MS0,u,beta,q5) {
i1   <- beta^2*(2*pnorm(u)-1)
i2   <- integrate(s2s2phi.n,lower=-u,upper=u)$value
i3   <- integrate(s2phi.n,lower=-u,upper=u)$value
Ex   <- as.matrix(xbar)
Q3   <- i1+i2+q5^2*MS0%*%E.S%*%t(MS0) #+q4^2*t(Ex)%*%MT0%*%E.S%*%t(MT0)%*%Ex
Q3   <- Q3 - 2*beta*i3 - 2*beta*q5*MS0%*%E0l #- 2*beta*q4*t(Ex)%*%MT0%*%E0l
Q3   <- Q3 + 2*q5*MS0%*%E2 #+ 2*q4*t(Ex)%*%MT0%*%E2 + 2*q4*q5*t(Ex)%*%MT0%*%E.S%*%t(MS0)
Q3}

