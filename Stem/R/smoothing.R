`smoothing` <-
function(ss) {
  m <- ss$m
  C <- ss$C
  m0<- ss$m0
  C0<- ss$C0
  d <- ss$d
  nobs <- ss$n

  mu <- matrix(NA, nobs, d)
  mu[nobs,] <- t(ss$Fmat) %*% m[nobs,]

  for (tt in (nobs-1):1)    {

	if (ss$p == 1)
		nextstep <-
			smootherstep.uni(
				m[tt,],
				C[[tt]],
				ss$Gmat,
				ss$Wmat,
				m[tt+1,],
				C[[tt+1]])

	else
		nextstep <-
			smootherstep(
				matrix(m[tt,],nrow=1),
				C[[tt]],
				ss$Gmat,
				ss$Wmat,
				matrix(m[tt+1,],nrow=1),
				C[[tt+1]])
                     
	m[tt,]  <- nextstep$ms
	C[[tt]] <- nextstep$Cs
	mu[tt,] <- t(ss$Fmat) %*% m[tt,]
	}
    
	ss$m0 <- m0
	ss$C0 <- C0
	if (is.ts(ss$Z))   ss$m <- ts(m,start(ss$z),end=end(ss$z),frequency=frequency(ss$z))
	else
	ss$m <- m
	ss$C <- C
	ss$mu <- mu
	ss
}

