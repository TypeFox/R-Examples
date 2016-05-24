grad <-
function(x,omega,l_period,n,freq,h){

 -colSums(kern(omega,h,freq)$v*as.vector((-1+exp(l_period-x[1]-t((t(freq)-omega))%*%x[2:3])))*cbind(rep(1,n),t((t(freq)-omega))))

}
