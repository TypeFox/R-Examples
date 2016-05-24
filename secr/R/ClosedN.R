# Various closed-population estimators
# Input is 'capthist' object of package 'secr'
# Murray Efford April 2010
# Based on Otis et al. 1978, C code of Anne Chao etc.
# modified 2010-09-05 for AICcwt
# modified 2011-05-04 for only non-ML estimators
# source('d:\\density secr 1.4\\secr\\r\\closedN.R')

############################################################################

closedN <- function (object, estimator = NULL, level = 0.95, maxN = 1e7, dmax = 10) {
    if (!inherits(object, 'capthist'))
        stop ("requires 'capthist' object")
    if (inherits(object, 'list'))
        lapply(object, closedN, estimator=estimator, level=level, maxN=maxN)
    else {
        allowed <- c('null','zippin','darroch','h2','betabinomial',
                     'jackknife','chao','chaomod','chao.th1', 'chao.th2')
        models <- c('M0','Mb','Mt','Mh','Mh', 'Mh','Mh','Mh','Mth','Mth')
        if (is.null(estimator)) estimator <- allowed
        else if ('ML' %in% estimator) estimator <- allowed[1:5]
        if (any(! estimator %in% allowed)) stop ("estimator not recognised")
        anyML <- any(estimator %in% allowed[1:5])
        ## Chao's Mth estimators need conventional capture histories
        if (any(estimator %in% c('chao.th1', 'chao.th2')))
            X <- apply(abs(object), 1:2, sum) > 0

        nocc <- ncol(object)
        counts <- as.matrix(summary(object)$counts)[,-(nocc+1)]
        getest <- function(x) {
            switch (x,
                null = null.est.c (counts['n',], nrow(object), maxN),
                zippin = zip.est.c (counts['n',], counts['u',], maxN),
                darroch = darr.est.c (counts['n',], nrow(object), maxN),
                jackknife = jackknife.est (counts['f',]),
                chao = chao.est (counts['f',]),
                chaomod = chaomod.est (counts['f',]),
                h2 = Mh2.est (counts['f',]),
                betabinomial = Mhbeta.est (counts['f',]),
                chao.th1 = chao.th1.est(X),
                chao.th2 = chao.th2.est(X)
            )
        }
        temp <- sapply (estimator,getest)
        tempCL <- apply(temp,2,LN.interval, level)
        temp <- as.data.frame(t(rbind(temp, tempCL)))
        temp[temp>maxN] <- NA
        novalue <- rep(NA, nrow(temp))

        temp$model <- models[match(estimator, allowed)]

        if (anyML) {
            temp$AIC <- -2*temp$LL + 2 * temp$npar
            temp$AICc <- ifelse ((temp$Mt1 - temp$npar - 1) > 0,
                -2 * temp$LL + 2 * temp$npar * temp$Mt1 / (temp$Mt1 - temp$npar - 1),
                NA)
            temp$dAICc <- temp$AICc-min(temp$AICc, na.rm=T)
            OK <- abs(temp$dAICc) < abs(dmax)
            sumdAICc <- sum(exp(-temp$dAICc[OK]/2), na.rm=T)
            temp$AICwt <- ifelse ( OK, round(exp(-temp$dAICc/2) / sumdAICc,4), 0)
        }
        else {
            ## 2011-05-04
            temp$AIC <- novalue
            temp$AICc <- novalue
            temp$dAICc <- novalue
            temp$AICwt <- novalue
        }

        ## added 2010-09-05

        temp <- temp[,c(8,4,5,9:12, 1:3,6,7)]
        names(temp) <- c('model','npar','loglik','AIC','AICc','dAICc', 'AICcwt', 'Mt1','Nhat',
                         'seNhat','lclNhat', 'uclNhat')
        temp[,3] <- round(as.matrix(temp[,3]), 3)
        temp[,c(4:6,9:12)] <- round(as.matrix(temp[,c(4:6,9:12)]), 2)
        temp[,7] <- round(as.matrix(temp[,7]), 3)
        temp
    }
}
############################################################################

#################################
## Null estimator M0
## continuous-N

null.est.c <- function (nt, Mt1, maxN=1e7) {
  loglik <- function (N) # Otis et al 1978 p105, modified for continuous N
                n*log(n) +
                (tt*N - n) * log (tt*N - n) -
                (tt*N) * log(tt*N) +
                lgamma (N+1) - lgamma(N-Mt1+1) - lgamma(Mt1+1)
  tt <- length(nt)
  n  <- sum(nt)
  fit <- optimize(f = loglik, interval = c(Mt1, maxN),  maximum = TRUE)
  nhat <- fit$maximum
  phat <- n / (tt*nhat)
  senhat <- sqrt(nhat / ( (1 - phat)^ -tt - tt/(1 - phat) + tt - 1) )
  c(Mt1 = Mt1, Nhat = nhat, seNhat = senhat, npar = 2, LL = fit$objective)
}

#################################
## Zippin estimator Mb
## continuous-N

zip.est.c <- function (nt, ut, maxN=1e7) {
   # Otis et al 1978 p108, modified for continuous N
   # mod 2010 04 03 for full likelihood
  loglik <- function (N) {
                 cterm <- ifelse (m>0, m*log(m/M.) + (M.-m) * log(1-m/M.), 0)
                 lgamma(N+1) - lgamma(N-Mt1+1) - lgamma(Mt1+1) +
                     Mt1*log(Mt1) +
                    (tt*N - M. - Mt1) * log(tt*N - M. - Mt1) -
                    (tt*N - M.) * log(tt*N - M.) +
                    cterm
            }
  tt <- length(ut)
  m <- sum(nt-ut)
  Mj <- c(0,cumsum(ut))
  M.  <- sum(Mj[2:tt])
  Mt1 <- Mj[tt+1]
  ## OK <- sum( (tt+1-2*(1:tt)) * ut) > 0
  fit <- optimize(f = loglik, interval = c(Mt1, maxN),  maximum = TRUE)
  nhat <- fit$maximum
  phat <- Mt1 / (tt*nhat - M.)
  senhat <- sqrt((nhat * (1-phat)^tt * (1 - (1-phat)^tt)) /
                 ((1 - (1-phat)^tt)^2 - tt^2*phat^2*(1-phat)^(tt-1)))
  c(Mt1 = Mt1, Nhat = nhat, seNhat = senhat, npar = 3, LL = fit$objective)
}

#################################
## Darroch estimator Mt
## continuous-N

darr.est.c <- function (nt, Mt1, maxN=1e7) {
  loglik <- function (N)  # Otis et al 1978 p106-7, modified for continuous N
                lgamma(N+1) - lgamma(N-Mt1+1) - lgamma(Mt1+1) +
                -tt * N * log(N) + sum (nt * log(nt) + (N-nt)*log(N-nt))
  tt <- length(nt)
  fit <- optimize(f = loglik, interval = c(Mt1, maxN),  maximum = TRUE)
  nhat <- fit$maximum
  phat <- sum(nt) / (tt*nhat)
  t1 <- prod(1-nt/nhat)
  t2 <- sum(1/(1-nt/nhat))
  senhat <- sqrt(nhat / (1/t1 + tt - 1 -t2))
  c(Mt1 = Mt1, Nhat = nhat, seNhat = senhat, npar = tt+1, LL = fit$objective)
}

#################################
## Chao Mh

chao.est <- function (fi) {
    tt <- length(fi)	# number of capture occasions
    S <- sum(fi)		# total number of individuals caught
    N <- S + fi[1]^2 / fi[2] / 2
    varN <- fi[2] * (0.25 * (fi[1]/fi[2])^4 + (fi[1]/fi[2])^3 + 0.5 * (fi[1]/fi[2])^2)
    c(Mt1 = S, Nhat = N, seNhat = varN^0.5, npar = NA, LL = NA)
}

#################################
## Chao modified Mh

chaomod.est <- function (fi) {
    tt <- length(fi)	# number of capture occasions
    S <- sum(fi)		# total number of individuals caught
    ## cat (tt*fi[1], 2*fi[2], tt*fi[2], 3*fi[3], 3*fi[1]*fi[3], 2*fi[2]^2, '\n')
    Nmod <- ifelse ((tt*fi[1] > 2*fi[2]) & (tt*fi[2] > 3*fi[3]) & (3*fi[1]*fi[3] > 2*fi[2]^2),
                    S + (fi[1]^2 / fi[2] / 2) * (1 - 2*fi[2]/
                    (tt*fi[1])) / (1 - 3 * fi[3] / (tt * fi[2])),NA)
    if (!is.na(Nmod)) {
	A	<- (1 - 2*fi[2]/(tt * fi[1])) / (1 - 3*fi[3]/(tt * fi[2]))
	varNmod	<- fi[2] * (0.25 * A^2 * (fi[1]/fi[2])^4 + A^2 *
                  (fi[1]/fi[2])^3 + 0.5 * A * (fi[1]/fi[2])^2)
    }
    else varNmod <- NA
    if (is.na(Nmod))
        chao.est(fi)
    else
        c(Mt1 = S, Nhat = Nmod, seNhat = varNmod^0.5, npar = NA, LL = NA)
}

#################################
## 2-part mixture

Mh2.est <- function (fi) {
  ## based in part on S+ code of Shirley Pledger 24/4/98
  loglik <- function (beta) {
    N <- Mt1 + exp(beta[1])
    pi1 <- invlogit (beta[2])
    theta1 <- invlogit (beta[3])
    theta2 <- invlogit (beta[4])
    i <- 1:tt
    terms <-  log(pi1    * theta1^i * (1-theta1)^(tt-i) +
                 (1-pi1) * theta2^i * (1-theta2)^(tt-i))
    LL <- lgamma (N+1) - lgamma(N-Mt1+1) - lgamma(Mt1+1) +
        (N-Mt1) * log (pi1 * (1-theta1)^tt + (1-pi1) * (1-theta2)^tt) +
         sum(fi*terms)
      -LL
  }
  tt <- length(fi)
  Mt1  <- sum(fi)
  start <- c(log(5), 0, -1, 1)
  fit <- nlm (p = start, f = loglik, hessian = TRUE)
  nhat <- exp(fit$estimate[1]) + Mt1
  pi1 <- invlogit(fit$estimate[2])
  theta1 <- invlogit(fit$estimate[3])
  theta2 <- invlogit(fit$estimate[4])
  mu <- pi1 * theta1 +  (1-pi1) * theta2
  vr <- pi1 * (theta1-mu)^2 + (1-pi1) * (theta2-mu)^2

    if (!is.null(fit$hessian)) {
        covar <- try(solve(fit$hessian))
        if (inherits(covar, "try-error")) {
            warning ("could not invert Hessian to compute ",
                     "variance-covariance matrix")
            covar <- matrix(rep(NA,16), ncol = 4)  # failed
        }
        senhat <-  exp(fit$estimate[1]) * sqrt(exp(covar[1,1]^2)-1)
    }
    else senhat <- NA
  c(Mt1 = Mt1, Nhat = nhat, seNhat = senhat, npar = 4, LL = -fit$minimum)
}

#########################################
## beta binomial
## Dorazio & Royle
## based in part on S+ code of Shirley Pledger 24/4/98

Mhbeta.est <- function (fi, maxN = 1e7) {
    loglik <- function (pr) {
        pr <- exp(pr)   ## all on log scale
        N <- Mt1 + pr[1]
        rat   <- pr[2] * (1- pr[2])/ pr[3]
        if ((rat<1) | (N>maxN))  return (1e10)
        alpha <- pr[2] * (rat-1)
        beta  <- (1 - pr[2]) * (rat-1)
        i <- 1:tt
        terms <-  lgamma(alpha+i) + lgamma(beta+tt-i)
        LL <- lgamma (N+1) - lgamma(N-Mt1+1) - lgamma(Mt1+1) +
            N * (lgamma(alpha+beta) - lgamma(alpha) - lgamma(beta) -
            lgamma(alpha + beta + tt)) +
           (N-Mt1) * (lgamma(alpha) + lgamma(beta+tt)) +
            sum(fi*terms)
        -LL
    }
    tt <- length(fi)
    Mt1  <- sum(fi)
    start <- log(c(10, 1/tt, 0.2 * 1/tt * (1 - 1/tt) ))
    fit <- nlm (p = start, f = loglik, hessian = TRUE)
    nhat <- exp(fit$estimate[1]) + Mt1

    if (!is.null(fit$hessian)) {
        covar <- try(solve(fit$hessian))
        if (inherits(covar, "try-error")) {
            warning ("could not invert Hessian to compute ",
                     "variance-covariance matrix")
            covar <- matrix(rep(NA,16), ncol = 4)  # failed
        }
        senhat <-  exp(fit$estimate[1]) * sqrt(exp(covar[1,1]^2)-1)
    }
    else senhat <- NA
  c(Mt1 = Mt1, Nhat = nhat, seNhat = senhat, npar = 3, LL = -fit$minimum)
}

#################################
## jackknife Mh

jackknife.est <- function (fi, full = F) {
# Calculate Burnham & Overton's jackknife estimate for closed populations
#
# inp is a vector of capture frequencies
#
# MODIFIED
# 30/3/95 Fix bug when mt=5 (last jacknife selected) - length(test) s/b 5
# 30/3/95 Optional full output
# 30/3/95 Add confIDence limits to 'short' output
# 3/4/95  Implement unconditional se using method of K.P.Burnham
# 3/4/95  Minor changes to full output
# 1/4/2010 Added test for illconditioned data N<S cf CAPTURE jack.for L351
# 1/4/2010 do not use deads
#

	jack.fill <- function(tt)
	{

	# Input:  tt is the number of capture occasions (e.g., days)
	# Output: matrix of jackknife coefficients for Burnham & Overton estimator

	# Murray Efford 30/3/95

		T1 <- tt - 1
		T2 <- tt - 2
		T3 <- tt - 3
		T4 <- tt - 4
		T5 <- tt - 5
		occ5 <- min(tt, 5)
		fcoeff <- matrix(data = 0, ncol = 5, nrow = 5)
		fcoeff[1, 1] <- T1/tt
		fcoeff[2, 1] <- (2 * tt - 3)/tt
		fcoeff[2, 2] <-  - T2^2/(tt * T1)
		fcoeff[3, 1] <- (3 * tt - 6)/tt
		fcoeff[3, 2] <-  - (3 * tt^2 - 15 * tt + 19)/(tt * T1)
		fcoeff[3, 3] <- T3^3/(tt * T1 * T2)
		fcoeff[4, 1] <- (4 * tt - 10)/tt
		fcoeff[4, 2] <-  - (6 * tt^2 - 36 * tt + 55)/(tt * T1)
		fcoeff[4, 3] <- (4 * tt^3 - 42 * tt^2 + 148 * tt - 175)/(tt * T1 * T2)
		fcoeff[4, 4] <-  - T4^4/(tt * T1 * T2 * T3)
		fcoeff[5, 1] <- (5 * tt - 15)/tt
		fcoeff[5, 2] <-  - (10 * tt^2 - 70 * tt + 125)/(tt * T1)
		fcoeff[5, 3] <- (10 * tt^3 - 120 * tt^2 + 485 * tt - 660)/(tt * T1 * T2)
		fcoeff[5, 4] <-  - (T4^5 - T5^5)/(tt * T1 * T2 * T3)
		fcoeff[5, 5] <- T5^5/(tt * T1 * T2 * T3 * T4)
		fcoeff <- fcoeff[1:occ5, 1:occ5]	# Use sub-matrix if tt < 5
		if(tt > 5) {
			fcoeff <- cbind(fcoeff, matrix(data = 0, nrow = 5, ncol = tt - 5))
		}
		fcoeff + 1	# Adding one allows: Nj<-fcoeff%*%fi
	}

        fi      <- unlist(fi)
	S	<- sum(fi)
	occ	<- length(fi)
	occ5	<- min(5, occ)
	aki	<- jack.fill(occ)
	nj		<- aki %*% fi
	varnj 	<- (aki^2 %*% fi) - nj
	difnk	<- nj[2:occ5] - nj[1:(occ5 - 1)]
	b2f	<- ((aki[2:occ5,  ] - aki[1:(occ5 - 1),  ])^2) %*% fi
	test	<- (difnk/sqrt(S/(S - 1) * (b2f - difnk^2/S)))^2
	test[is.na(test)] <- 0
	test	<- rbind(test, 0)
	pk		<- 1 - pchisq((difnk/sqrt(S/(S - 1) * (b2f - difnk^2/S)))^2, 1)
	pk[difnk == 0] <- 1

	# Select jacknife on basis of test results, and interpolate estimate
	# The conditional variance calculations commented out here are superceded by the
	# unconditional calculations below

            first	<- function(vec) match(1, vec)
	     mt		<- first(test < 3.84) #
	     if(mt == 1) {
	          xtest	<- (nj[1] * fi[1])/(nj[1] - fi[1])
	          if(xtest > 3.84) {
	               alpha	<- (xtest - 3.84)/(xtest - test[1])
	               beta	<- 1 - alpha
	               z	<- aki[1,  ] * alpha + beta
	               N <- z %*% fi  # varN <- (z * z) %*% fi - N
	          }
	          else {
	               N	<- nj[1]
	               varN	<- varnj[1]
	          }
	     }
	     else {
	          alpha	<- (test[mt - 1] - 3.84)/(test[mt - 1] - test[mt])
	          beta	<- 1 - alpha
	          z 	<- aki[mt,  ] * alpha + aki[mt - 1,  ] * beta
	          N 	<- z %*% fi  # varN <- ((z * z) %*% fi) - N
	     }
             if (N<S) {
                 N <- nj[1]
                 varN <- NA
                 Pi <- NA
                 warning ("data illconditioned for jackknife")
             }
             else {
                 k1 <- occ5 - 1
        	 Pi <- rep(1, occ5)
        	 Beta <- pchisq(rep(3.8415, occ5), 1, pmax(test - 1, 0))
        	 for(i in 2:occ5) Pi[i] <- Pi[i - 1] * (1 - Beta[i - 1])
        	 Pi[1:k1] <- Pi[1:k1] * Beta[1:k1]
        	 if(occ5 < 5) Pi[occ5] <- 1 - sum(Pi[1:k1])
                 varN <- sum(Pi * varnj) + sum(Pi * (nj - sum(Pi * nj))^2)
              }

	# Output
             estimates <-  c(Mt1 = sum(fi), Nhat = N, seNhat = sqrt(varN), npar = NA, LL = NA)
	     if(full) list(fi = fi, estimates = estimates, mt = mt, aki = aki,
	               nj = data.frame(nj, se = sqrt(varnj), X2 = test, Pi = Pi, Beta = Beta))
	     else estimates
}
############################################################################

chao.th1.est <- function (X) {
    ## Translated from Anne Chao's C code via Pascal by MGE
    ## Code from Anne Chao 17/3/03
    ## X is a matrix representing capture histories
    tt <- ncol(X)
    D <- nrow(X)
    ni <- apply(X,2,sum)
    Y <- apply(X,1,sum)
    fi <- tabulate (Y, nbins = tt)
    i <- 1:tt
    temp1 <- sum (i * fi)
    temp2 <- sum (i * (i-1)*fi)
    temp3 <- 0;
    for (i in 1:(tt-1))
        for (j in (i+1):tt)
            temp3 <- temp3 + 2*ni[i]*ni[j]
    sq2 <- temp1^-2
    c.1 <- 1 - fi[1]/temp1
    if (c.1 > 0) {
        n0  <- D/c.1     # fails if no recaptures
        r1  <- max(n0*temp2/temp3-1, 0.0)
        n1  <- n0 + fi[1]/c.1*r1
    }
    else n1 <- NA    ## cannot compute

    v_th1 <- NA
    if (! is.na(n1)) {
        dn0 <- numeric(tt)
        dn0[1] <- 1/c.1 * (1 + D/temp1)
        dn0[2:tt] <- 1/c.1 - c.1^-2 * (2:tt) * fi[1]*sq2*D
        if (abs(r1) < 1e-6) {
            sum1 <- sum( dn0^2 * fi)
            sum2 <- sum( dn0 * fi)
            v_th1 <- (sum1-sum2^2/n1);
        }
        else {
            dr1 <- dn0[1] * temp2 / temp3
            dnf <- numeric(tt)
            dnf[1] <- dn0[1] + r1/c.1 + fi[1] * r1 / c.1 / temp1 + fi[1] / c.1 * dr1
            j <- 2:tt
            drj <- dn0[j] * temp2 / temp3 + n0 * j * (j-1) / temp3
            dnf[j] <- dn0[j] - fi[1] * r1 * c.1^-2 * j * fi[1] * sq2 + fi[1] / c.1 * drj
            dnk <- fi[1] / c.1 * (-2 * n0 * temp2 * temp3^-2 * (temp1-ni))
            sum1 <- sum (dnf^2 * fi)
            sum2 <- sum (dnf * fi)
            su1 <- sum( dnk^2 * ni)
            su2 <- sum( dnk * ni)
            su3 <- 0;

            nk_nl <- function (k,l) {
                sum (X[1:D,k]==1 & X[1:D,l]==1)
            }
            for (i in 1:(tt-1))
                for (j in (i+1):tt)
                    su3 <- su3 + 2*dnk[i]*dnk[j]*nk_nl(i,j)
            fk_nl <- function (k,l) {
                sum ((Y[1:D]==k) * X[1:D,l])
            }
            su4 <- 0;
            for (i in 1:tt)
                for (j in 1:tt)
                    su4 <- su4 + 2*dnf[i]*dnk[j]*(fk_nl(i,j)- fi[i] *  ni[j] / n1)
            v_th1 <- (sum1 - sum2 * sum2 / n1) + (su1 - su2 * su2 / n1) + su3 + su4
        }
    }
    c(Mt1 = D, Nhat = n1, seNhat = v_th1^0.5, npar = NA, LL = NA)
}
############################################################################

chao.th2.est <- function (X) {
    ## Translated from Anne Chao's C code via Pascal by MGE
    ## Code from Anne Chao 17/3/03
    ## X is a matrix representing capture histories
    tt <- ncol(X)
    D <- nrow(X)
    ni <- apply(X,2,sum)
    Y <- apply(X,1,sum)
    fi <- tabulate (Y, nbins = tt)
    i <- 1:tt
    temp1 <- sum (i * fi)
    temp2 <- sum (i * (i-1)*fi)
    temp3 <- 0;
    for (i in 1:(tt-1))
        for (j in (i+1):tt)
            temp3 <- temp3 + 2*ni[i]*ni[j]
    sq2 <- temp1^-2

    c1 <- fi[1] - 2*fi[2]/(tt-1)
    c.2 <- 1 - c1/temp1
    if (c.2 > 1e-6){
    ## should probably restrict as follows 2013-04-20
    ## if ((c.2 > 1e-6) & (tt>3)){
        n0  <- D/c.2     # fails if no recaptures
        r2  <- max(n0*temp2/temp3-1, 0.0)
        n2  <- n0 + fi[1]/c.2*r2
    }
    else n2 <- NA    ## cannot compute

    v_th2 <- NA
    if (! is.na(n2)) {
        dn0 <- numeric(tt)
        dn0[1] <- 1/c.2 * (1 + D/temp1)
        dn0[2] <- 1/c.2 - 2*c.2^-2 * (1/temp1/(tt-1)+c1*sq2)*D
        if (tt>2)  ## inserted 2013-04-20 to eliminate warning
        dn0[3:tt] <- 1/c.2 - c.2^-2 * (3:tt) * c1 *sq2 * D
        if (abs(r2) < 1e-6) {
            sum1 <- sum( dn0^2 * fi)
            sum2 <- sum( dn0 * fi)
            v_th2 <- (sum1-sum2^2/n2)
        }
        else {
            dr1 <- dn0[1] * temp2 / temp3
            dnf <- numeric(tt)
            dnf[1] <- dn0[1] + r2/c.2 + fi[1] * r2 / c.2 / temp1 + fi[1] / c.2 * dr1
            dr2    <- dn0[2] * temp2 / temp3 + 2 * n0 / temp3
            dnf[2] <- dn0[2] - 2 * fi[1] * r2 * c.2^-2 * ( 1/temp1/(tt-1) + c1 * sq2) +
                fi[1] / c.2 * dr2
            j <- 3:tt
            drj    <- dn0[j] * temp2/temp3 + n0 * j * (j-1)/temp3
            dnf[j] <- dn0[j] - fi[1] * r2* c.2^-2 * j * c1 * sq2 + fi[1] / c.2 * drj
            dnk <- fi[1] / c.2 * (-2 * n0 * temp2 * temp3^-2 * (temp1-ni))
            sum1 <- sum (dnf^2 * fi)
            sum2 <- sum (dnf * fi)
            su1 <- sum( dnk^2 * ni)
            su2 <- sum( dnk * ni)
            su3 <- 0;

            nk_nl <- function (k,l) {
                sum (X[1:D,k]==1 & X[1:D,l]==1)
            }
            for (i in 1:(tt-1))
                for (j in (i+1):tt)
                    su3 <- su3 + 2*dnk[i]*dnk[j]*nk_nl(i,j)
            fk_nl <- function (k,l) {
                sum ((Y[1:D]==k) * X[1:D,l])
            }
            su4 <- 0;
            for (i in 1:tt)
                for (j in 1:tt)
                    su4 <- su4 + 2*dnf[i]*dnk[j]*(fk_nl(i,j)- fi[i] *  ni[j] / n2)
            v_th2 <- (sum1 - sum2 * sum2 / n2) + (su1 - su2 * su2 / n2) + su3 + su4
        }
    }
    c(Mt1 = D, Nhat = n2, seNhat = v_th2^0.5, npar = NA, LL = NA)
}
############################################################################

LN.interval <- function (x, level) {
     # Confidence limits as in Rexstad & Burnham 1991 CAPTURE Users' Guide
     Mt1 <- x[1]; Nhat <- x[2]; seNhat <- x[3]
     f0  <- Nhat - Mt1
     z <- qnorm (1 - (1-level)/2)
     cc1 <- ifelse (f0>0, exp(z * sqrt(log(1 + (seNhat/f0)^2))), NA)
     c(lcl = Mt1 + f0 / cc1, ucl = Mt1 + f0 * cc1)
}

############################################################################

