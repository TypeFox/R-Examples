hess <-
function(x,omega,l_period,n,freq,h){

 he=colSums(kern(omega,h,freq)$v*as.vector(exp(l_period-x[1]-t((t(freq)-omega))%*%x[2:3]))*cbind(rep(1,n),t((t(freq)-omega)),t((t(freq)-omega))[,1]*t((t(freq)-omega))[,2],t((t(freq)-omega))^2))
 matrix(c(he[1:3],he[c(2,5,4)],he[c(3,4,6)]),3,3)

}
