`filtering` <-
function(ss) {
	m <- matrix(NA,ss$n,ss$p)
	C <- vector("list",ss$n)
	firststep <-
		filterstep(
			z = matrix(ss$z[1,]),
			Fmat	= ss$Fmat,
			Gmat = ss$Gmat,
			Vt 	= ss$Vmat,
			Wt	= ss$Wmat,
			mx	= ss$m0,
			Cx	= ss$C0,
			XXXcov  = ss$XXX[,,1],
			betacov = ss$beta,
			flag    = ss$flag.cov
               )
	m[1,] <- firststep$m
	C[[1]]<- firststep$C
	loglik<- firststep$loglikterm

  ## run the recursion
	for (tt in 2:ss$n) {
	nextstep <-
		filterstep(
			z= matrix(ss$z[tt,]),
			Fmat = ss$Fmat,
			Gmat = ss$Gmat,
			Vt = ss$Vmat,
			Wt = ss$Wmat,
			mx 	= matrix(m[tt-1,],nrow=1),
			Cx 	= C[[tt-1]],
			XXXcov  = ss$XXX[,,tt],
			betacov = ss$beta,
			flag    = ss$flag.cov
		)
      m[tt,]  <- nextstep$m
      C[[tt]] <- nextstep$C
      loglik  <- loglik + nextstep$loglikterm
   }
   
	if (is.ts(ss$z))
		ss$m <- ts(m,start(ss$z),end=end(ss$z),frequency=frequency(ss$z))
	else
	ss$m <- m
	ss$C <- C
	ss$loglik <- loglik
	
	ss
}

