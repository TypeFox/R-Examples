lagrange <-
function(gamm0, nstrata , nn1 , nn0, n1a , n0a , grpa ,repp , xx , ofs , yy)
{
	fv <- xx%*%gamm0 + ofs 
	fv <- exp(fv)
	fv <- fv/(1+fv)
	r  <- fv[1:nstrata]
	r  <- nn1 - (nn1 + nn0)*r
	r  <- c(0,r)
	r1 <- (n1a-r)*n0a
	r2 <- (n0a+n1a)*r
	g  <- r1[grpa]*fv
	g  <- g/(r1[grpa] + r2[grpa]*(1-fv))
	g  <- t(xx)%*%((1-g)*yy - g*(repp-yy)) 
	g
}
