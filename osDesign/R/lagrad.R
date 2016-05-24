lagrad <-
function(gamm0, nstrata, nobs , ncovs ,nn1 , nn0 , n1a , n0a ,ee ,repp , grpa , n0 , n1 , xx , ofs)
{
	fv  <- xx%*%gamm0 + ofs
	fv  <- exp(fv)
	fv  <- fv/(1+fv)
	xxx <- xx[(nstrata+1):(nstrata+nobs),]
	grp <- grpa[(nstrata+1):(nstrata+nobs)]-1
	rep <- repp[(nstrata+1):(nstrata+nobs)]
	e   <- ee[(nstrata+1):(nstrata+nobs),2:(nstrata+1)]
	pp  <- fv[1:nstrata]
	g   <- fv[(nstrata+1):(nstrata+nobs)]
	r   <- nn1 - (nn1+nn0)*pp
	r0  <- r*(n0+n1)
	r1  <- n0*n1*(nn0+nn1)*(n0+n1)*pp*(1-pp)
	r2  <- n0*(n1-r)

	d <- r1[grp]*g*(1-g)*rep/(r2[grp]+r0[grp]*(1-g))^2
  d <- as.vector(d)                                       ## ADDITIONAL 29TH OCT 2007

	de <- d*e
	g  <- -t(xxx)%*%de
	g  <- cbind(g,matrix(0,(nstrata+ncovs),ncovs))
	r  <- c(0,r)
	r1 <- n0a*(n1a-r)
	r2 <- (n0a+n1a)*r
	r0 <- n0a*n1a*(n1a-r)*(n0a+r)
	fv <- r0[grpa]*fv*(1-fv)*repp/(r1[grpa]+r2[grpa]*(1-fv))^2
	fv <- t(xx)%*%(hdp(fv,xx))
	g <- g - fv
	g
}
