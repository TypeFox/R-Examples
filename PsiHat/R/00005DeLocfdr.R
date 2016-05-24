##from the R package locfdr (v1.1.7) by Bradley Efron, Brit B. Turnbull, and Balasubramanian Narasimhan.

locmle <- function(z, xlim , Jmle = 35, d = 0, s = 1, ep = 1/100000, sw = 0, Cov.in)
{
        ## uses z-values in [-xlim,xlim] to find mles for p0,del0,sig0
        ## Jmle number of iterations, beginning at (del0,sig0)=(d,s)
        ## sw=1 returns correlation matrix
        N = length(z)
        if (missing(xlim)) {
          if (N>500000) b = 1
          else b=4.3 * exp(-0.26*log(N,10)) 
          xlim=c(median(z),b*diff(quantile(z)[c(2,4)])/(2*qnorm(.75)))
        }
        aorig=xlim[1]-xlim[2]
        borig=xlim[1]+xlim[2]
        z0=z[which(z>=aorig & z<=borig)]
        N0 = length(z0)
        Y = c(mean(z0), mean(z0^2))
        that = N0/N
        ######## find mle estimates ###########################
        for(j in 1:Jmle) {
                bet = c(d/s^2, -1/(2 * s^2))
                aa = (aorig - d)/s
                bb = (borig - d)/s
                H0 = pnorm(bb) - pnorm(aa)
                fa = dnorm(aa)
                fb = dnorm(bb)
                H1 = fa - fb
                H2 = H0 + aa * fa - bb * fb
                H3 = (2 + aa^2) * fa - (2 + bb^2) * fb
                H4 = 3 * H0 + (3 * aa + aa^3) * fa - (3 * bb + bb^3) * fb
                H = c(H0, H1, H2, H3, H4)
                r = d/s
                I = matrix(rep(0, 25), 5)
                for(i in 0:4)
                        I[i + 1, 0:(i + 1)] = choose(i, 0:i)
                u1 = s^(0:4)
                II = pmax(row(I) - col(I), 0)
                II = r^II
                I = u1 * (I * II)
                E = as.vector(I %*% H)/H0
                E1 = E[2]
                E2 = E[3]
                E3 = E[4]
                E4 = E[5]
                mu = c(E1, E2)
                V = matrix(c(E2 - E1^2, E3 - E1 * E2, E3 - E1 * E2, E4 -
                        E2^2), 2)
                bett = bet + solve(V, Y - mu)/(1 + 1/j^2)
                if(bett[2]>0)bett=bet+.1*solve(V, Y - mu)/(1 + 1/j^2)
                if (is.na(bett[2])) break
                else if (bett[2]>=0) break
                d =  - bett[1]/(2 * bett[2])
                s = 1/sqrt(-2 * bett[2])
                if(sum((bett - bet)^2)^0.5 < ep)
                        break
        }
        if (is.na(bett[2])) {
          mle = rep(NA, 6)
          Cov.lfdr = NA
          Cor = matrix(NA, 3, 3)
        }
        else if (bett[2] >=0)  {
          mle = rep(NA, 6)
          Cov.lfdr = Cov.out = NA
          Cor = matrix(NA, 3, 3)
        }
        else {
          aa = (aorig - d)/s
          bb = (borig - d)/s
          H0 = pnorm(bb) - pnorm(aa)
          p0 = that/H0
          ############  sd calcs ###########################
          J = s^2 * matrix(c(1, 0, 2 * d, s), 2)
          JV = J %*% solve(V)
          JVJ = JV %*% t(J)
          mat2 = cbind(0, JVJ/N0)
          mat1 = c((p0 * H0 * (1 - p0 * H0))/N, 0, 0)
          mat = rbind(mat1, mat2)
          h = c(H1/H0, (H2 - H0)/H0)
          matt = c(1/H0,  - (p0/s) * t(h))
          matt = rbind(matt, cbind(0, diag(2)))
          C = matt %*% (mat %*% t(matt))
          mle = c(p0, d, s, diag(C)^0.5)
          if(sw == 1) {
                sd = mle[4:6]
                Co = C/outer(sd, sd)
                dimnames(Co) = list(c("p0", "d", "s"), c("p0", "d", "s"))
                Cor = Co[c(2, 3, 1), c(2, 3, 1)]
          }
          if (!missing(Cov.in)) {
            i0 = which(Cov.in$x > aa & Cov.in$x < bb)
            Cov.out = loccov(N, N0, p0, d, s, Cov.in$x, Cov.in$X, Cov.in$f,
                             JV, Y, i0, H, h, Cov.in$sw)
          }
        }
        names(mle) = c("p0", "del0", "sig0", "sd.p0", "sd.del0", "sd.sig0")
        mle = mle[c(2, 3, 1, 5, 6, 4)]
        out = list(mle=mle)
        if (sw==1) {
          Cor = list(Cor = Cor)
          out = c(out, Cor)
        }
        if (!missing(Cov.in)) {
          if (Cov.in$sw == 2) {
            pds. = list(pds.=Cov.out)
            out = c(out, pds.)
          }
          else if (Cov.in$sw == 3) {
            Ilfdr = list(Ilfdr=Cov.out)
            out = c(out, Ilfdr)
          }
          else {
            Cov.lfdr = list(Cov.lfdr=Cov.out)
            out = c(out, Cov.lfdr)
          }
        }
        if ((sw==1) | !missing(Cov.in)) return(out)
        else return(mle)
}

loccov <- function(N, N0, p0, d, s, x, X, f, JV, Y, i0, H, h, sw) {
  M = rbind(1, x - Y[1], x^2 - Y[2])
  if (sw==2) {
    K = length(x)
    K0 = length(i0)    
    toprow = c(1 - N0/N, -t(h) %*% JV / s)
    botrow = cbind(0, JV / p0)
    mat = rbind(toprow, botrow)
    M0 = M[,i0]
    dpds.dy0 = mat %*% M0 / N / H[1]
    dy0.dy = matrix(0, K0, K)
    dy0.dy[,i0] = diag(1, K0)
    dpds.dy = dpds.dy0 %*% dy0.dy
    rownames(dpds.dy) = c("p", "d", "s")
    return(dpds.dy)
  }
  else {
    xstd = (x - d)/s
    U = cbind(xstd - H[2]/H[1], xstd^2 - H[3]/H[1])
    M[,-i0] = 0
    dl0plus.dy = cbind(1 - N0/N, U %*% JV / s) %*% M /N/H[1]/p0  
    G <- t(X) %*% (f * X)
    dl.dy = X %*% solve(G) %*% t(X)
    dlfdr.dy = dl0plus.dy - dl.dy
    if (sw==3) return(dlfdr.dy)
    else {
      Cov.lfdr = dlfdr.dy %*% (f * t(dlfdr.dy))
      return(Cov.lfdr)
    }
  }
}

loccov2 <- function(X, X0, i0, f, ests, N) {
        d = ests[1]
        s = ests[2]
        p0 = ests[3]
        theo = I(ncol(X0)==1)
        Xtil <- X[i0,]
        X0til <- X0[i0,]
        G <- t(X) %*% (f * X)
        G0 <- t(X0til) %*% X0til
        B0 <- X0 %*% (solve(G0) %*% t(X0til)) %*% Xtil
        C <- B0 - X
        Ilfdr = C %*% solve(G, t(X))
        Cov <- C %*% solve(G) %*% t(C)
        if (theo)
          D = matrix(1,1,1)
        else
          D = matrix(c(1, 0, 0, d, s^2, 0, s^2 + d^2, 2 * d * s^2, s^3), 3)
	gam. = solve(G0, t(X0til)) %*% (Xtil %*% solve(G, t(X)))
	pds. = D %*% gam.
        if (theo) pds. = rbind(pds., matrix(0, 2, nrow(X)))
	pds.[1,] = pds.[1,] - 1/N
	m1 = pds. %*% f
	m2 = pds.^2 %*% f
	stdev = as.vector(sqrt(m2 - m1^2/N))
	stdev[1] = p0 * stdev[1]
        pds.[1,] = p0 * pds.[1,]
        rownames(pds.) = c("p", "d", "s")
        list(Ilfdr = Ilfdr, pds.=pds., stdev=stdev, Cov=Cov)
}



locfdr <-function(zz, bre = 120, df = 7, pct = 0, pct0 = 1/4, nulltype = 1, type = 0, plot = 1, mult, mlests, main = " ", sw = 0)
{
  call = match.call()
	if(length(bre) > 1) {
		lo <- min(bre)
		up <- max(bre)
		bre <- length(bre)
	}
	else {
		if(length(pct) > 1) {
			lo <- pct[1]
			up <- pct[2]
		}
		else {
			if(pct == 0) {
				lo <- min(zz)
				up <- max(zz)
			}
			if(pct < 0) {
				med = median(zz)
				ra = med + (1 - pct) * (range(zz) - med)
				lo = ra[1]
				up = ra[2]
			}
			if(pct > 0) {
				v <- quantile(zz, c(pct, 1 - pct))
				lo <- v[1]
				up <- v[2]
			}
		}
	}
	zzz <- pmax(pmin(zz, up), lo)
	breaks <- seq(lo, up, length = bre)
	zh <- hist(zzz, breaks = breaks, plot = F)
	x <- (breaks[-1] + breaks[ - length(breaks)])/2
	yall <- y <- zh$counts
	K <- length(y)
	N <- length(zz)
	if(pct > 0) {
		y[1.] <- min(y[1.], 1.)
		y[K] <- min(y[K], 1.)
	}
	if(type == 0) {
                X <- cbind(1, ns(x, df = df))
		f <- glm(y ~ ns(x, df = df), poisson)$fit
	}
	else {
                X <- cbind(1, poly(x, df = df))
		f <- glm(y ~ poly(x, df = df), poisson)$fit
	}
	l <- log(f)
	Fl <- cumsum(f)
	Fr <- cumsum(rev(f))
	D <- (y - f)/(f + 1)^0.5
	D <- sum(D[2:(K - 1)]^2)/(K - 2 - df)
	if(D > 1.5)
          warning(paste("f(z) misfit = ", round(
			D, 1), ".  Rerun with increased df", sep=""))
        # ............. create fp0 matrix ..........................
        if (nulltype == 3) {
                fp0 = matrix(NA, 6, 4)
                colnames(fp0) = c("delta", "sigleft", "p0", "sigright")
        }
        else {
                fp0 = matrix(NA, 6, 3)
                colnames(fp0) = c("delta", "sigma", "p0")
        }
        rownames(fp0) = c("thest", "theSD", "mlest", "mleSD", "cmest", "cmeSD")
        fp0["thest", 1:2] = c(0,1)
        fp0["theSD", 1:2] = 0
	# ..............begin central matching f0 calcs...............        
	imax <- seq(l)[l == max(l)][1]
	xmax <- x[imax]
	if(length(pct0) == 1) {
		pctup <- 1 - pct0
		pctlo <- pct0
	}
	else {
		pctlo <- pct0[1]
		pctup <- pct0[2]
	}
	lo0 <- quantile(zz, pctlo)
	hi0 <- quantile(zz, pctup)
	nx <- length(x)
	i0 <- (1.:nx)[x > lo0 & x < hi0]
	x0 <- x[i0]
	y0 <- l[i0]
	if(nulltype == 3) {
		X00 <- cbind((x0 - xmax)^2, pmax(x0 - xmax, 0)^2)
	}
	else {
		X00 <- cbind(x0 - xmax, (x0 - xmax)^2)
	}
	lr <- lm(y0 ~ X00)
	co <- lr$coef
        ## Error messages for failed CM estimation ##
        if (nulltype == 3) {
          cmerror = I(is.na(co[3]) | is.na(co[2]))
          if (!cmerror) cmerror = I(co[2] >= 0 | co[2]+co[3]>=0)
        }
        else {
          cmerror = is.na(co[3])
          if (!cmerror) cmerror = I(co[3] >= 0)
        }
        if (cmerror) {
          if (nulltype == 3)
            stop("CM estimation failed.  Rerun with nulltype = 1 or 2.")
          else
            if (nulltype == 2)
            stop("CM estimation failed.  Rerun with nulltype = 1.")
          else {
            X0 <- cbind(1, x - xmax, (x - xmax)^2)
            warning("CM estimation failed, middle of histogram non-normal")
          }
        }
        else {
	  if(nulltype == 3) {
		X0 <- cbind(1, (x - xmax)^2, pmax(x - xmax, 0)^2)
                sigs <- 1/sqrt(-2 * (c(co[2], co[2] + co[3])))
                fp0["cmest", c(1,2,4)] <- c(xmax, sigs)
	  }
	  else {
		X0 <- cbind(1, x - xmax, (x - xmax)^2)
                xmaxx <-  - co[2.]/(2. * co[3.]) + xmax
                sighat <- 1./sqrt(-2. * co[3.])
                fp0["cmest", 1:2] <- c(xmaxx, sighat)
	  }
	  l0 <- as.vector(X0 %*% co)
	  f0 <- exp(l0)
	  p0 <- sum(f0)/sum(f)
	  f0 <- f0/p0
          fp0["cmest", 3] <- p0
        }
	#............... begin MLE f0 calcs ........................
        b = 4.3 * exp(-0.26*log(N,10))
        if(missing(mlests)){
          med = median(zz);sc=diff(quantile(zz)[c(2,4)])/(2*qnorm(.75))
          mlests = locmle(zz, xlim=c(med, b*sc))
          if (N>500000) {
            warning("length(zz) > 500,000: For ML estimation, a wider interval than optimal was used.  To use the optimal interval, rerun with mlests = c(", mlests[1], ", ", b * mlests[2], ").\n", sep="")
            mlests = locmle(zz, xlim=c(med, sc))
          }
        }
	if (!is.na(mlests[1])) {
          if (N>500000) b = 1
          if (nulltype == 1) {
              Cov.in = list(x=x, X=X, f=f, sw=sw)
              ml.out = locmle(zz, xlim = c(mlests[1], b * mlests[2]),
                d=mlests[1], s=mlests[2], Cov.in=Cov.in)
              mlests = ml.out$mle
            }
            else  mlests = locmle(zz, xlim = c(mlests[1], b * mlests[2]),
                d=mlests[1], s=mlests[2])
            fp0["mlest", 1:3] = mlests[1:3]
            fp0["mleSD", 1:3] = mlests[4:6]
        }
        if (sum(is.na(fp0[c(3,5),1:2])) == 0 & nulltype > 1)
          if(abs(fp0["cmest",1] - mlests[1]) > 0.050000000000000003 |
             abs(log(fp0["cmest",2]/mlests[2])) > 0.050000000000000003)
		warning("Discrepancy between central matching and maximum likelihood estimates.\nConsider rerunning with nulltype = 1")
        ## Error messages for failed ML estimation ##
        if (is.na(mlests[1])) {
          if (nulltype == 1) {
            if (is.na(fp0["cmest", 1]))
              stop("CM and ML Estimation failed, middle of histogram non-normal")
            else stop("ML estimation failed.  Rerun with nulltype=2")
          }
          else warning("ML Estimation failed")
        }
	if(nulltype < 2) {
		delhat = xmax = xmaxx = mlests[1]
		sighat = mlests[2]
		p0 = mlests[3]
		f0 = dnorm(x, delhat, sighat)
		f0 = (sum(f) * f0)/sum(f0)
	}
        fdr = pmin((p0 * f0)/f, 1)
	f00 <- exp( - x^2/2)
	f00 <- (f00 * sum(f))/sum(f00)
	p0theo <- sum(f[i0])/sum(f00[i0])
        fp0["thest", 3] = p0theo
	fdr0 <- pmin((p0theo * f00)/f, 1)
	f0p <- p0 * f0
	if(nulltype == 0)
		f0p <- p0theo * f00
	F0l <- cumsum(f0p)
	F0r <- cumsum(rev(f0p))
	Fdrl <- F0l/Fl
	Fdrr <- rev(F0r/Fr)
	Int <- (1 - fdr) * f * (fdr < 0.90000000000000002)
	##### raise fdr to 1 near xmax .............
        if (sum(x <= xmax & fdr == 1) > 0)
          xxlo <- min(x[x <= xmax & fdr == 1])
        else xxlo = xmax
        if (sum(x >= xmax & fdr == 1) > 0)
          xxhi <- max(x[x >= xmax & fdr == 1])
        else xxhi = xmax
        if (sum(x >= xxlo & x <= xxhi) > 0)
          fdr[x >= xxlo & x <= xxhi] <- 1
        if (sum(x <= xmax & fdr0 == 1) > 0)
          xxlo <- min(x[x <= xmax & fdr0 == 1])
        else xxlo = xmax
        if (sum(x >= xmax & fdr0 == 1) > 0)
          xxhi <- max(x[x >= xmax & fdr0 == 1])
        else xxhi = xmax
	if (sum(x >= xxlo & x <= xxhi) > 0)
          fdr0[x >= xxlo & x <= xxhi] <- 1
	##################### raise fdr to 1 for mle option
	if(nulltype == 1) {
		fdr[x >= mlests[1] - mlests[2] & x <= mlests[1] + mlests[
			2]] = 1
		fdr0[x >= mlests[1] - mlests[2] & x <= mlests[1] + mlests[
			2]] = 1
	}
	p1 <- sum((1 - fdr) * f)/N
	p1theo <- sum((1 - fdr0) * f)/N
	fall <- f + (yall - y)
	########Efdr1 calculations
	Efdr <- sum((1 - fdr) * fdr * fall)/sum((1 - fdr) * fall)
	Efdrtheo <- sum((1 - fdr0) * fdr0 * fall)/sum((1 - fdr0) * fall)
	iup <- (1:K)[x >= xmax]
	ido <- (1:K)[x <= xmax]
	Eleft <- sum((1 - fdr[ido]) * fdr[ido] * fall[ido])/sum((1 - fdr[
		ido]) * fall[ido])
	Eleft0 <- sum((1 - fdr0[ido]) * fdr0[ido] * fall[ido])/sum((1 -
		fdr0[ido]) * fall[ido])
	Eright <- sum((1 - fdr[iup]) * fdr[iup] * fall[iup])/sum((1 - fdr[
		iup]) * fall[iup])
	Eright0 <- sum((1 - fdr0[iup]) * fdr0[iup] * fall[iup])/sum((
		1 - fdr0[iup]) * fall[iup])
	Efdr <- c(Efdr, Eleft, Eright, Efdrtheo, Eleft0, Eright0)
        Efdr[which(is.na(Efdr))] = 1
	names(Efdr) <- c("Efdr", "Eleft", "Eright", "Efdrtheo", "Eleft0",
		"Eright0")
	if(nulltype == 0)
		f1 <- (1 - fdr0) * fall
	else f1 <- (1 - fdr) * fall
	############ multiple sample size Efdr1 calculation
	if(!missing(mult)) {
		mul = c(1, mult)
		EE = rep(0, length(mul))
		for(m in 1:length(EE)) {
			xe = sqrt(mul[m]) * x
			f1e = approx(xe, f1, x, rule = 2, ties=mean)$y
			f1e = (f1e * sum(f1))/sum(f1e)
			f0e = f0
			p0e = p0
			if(nulltype == 0) {
				f0e = f00
				p0e = p0theo
			}
			fdre = (p0e * f0e)/(p0e * f0e + f1e)
			EE[m] = sum(f1e * fdre)/sum(f1e)
		}
		EE = EE/EE[1]
		names(EE) = mul
	}
	#................. Accuracy Calcs .................................
        Cov2.out = loccov2(X, X0, i0, f, fp0["cmest",], N)
        Cov0.out = loccov2(X, matrix(1, length(x), 1), i0, f, fp0["thest",], N)
        if(sw == 3) {
                if (nulltype==0) Ilfdr = Cov0.out$Ilfdr
                else if (nulltype==1) Ilfdr = ml.out$Ilfdr
		else if (nulltype==2) Ilfdr = Cov2.out$Ilfdr
                else stop("With sw=3, nulltype must equal 0, 1, or 2.")
	        return(Ilfdr)
              }
        if (nulltype == 0) Cov = Cov0.out$Cov
        else if (nulltype == 1) Cov = ml.out$Cov.lfdr
        else Cov = Cov2.out$Cov
	lfdrse <- diag(Cov)^0.5
        fp0["cmeSD",1:3] = Cov2.out$stdev[c(2,3,1)]
        if (nulltype==3) fp0["cmeSD",4] = fp0["cmeSD",2]
        fp0["theSD",3] = Cov0.out$stdev[1]
	########### sw==2 returns Influence function pds. ##########
	if(sw == 2) {
          if (nulltype==0) {
             pds = fp0["thest", c(3,1,2)]
             stdev = fp0["theSD", c(3,1,2)]
             pds. = t(Cov0.out$pds.)
          }
          else if (nulltype==1) {
            pds = fp0["mlest",c(3,1,2)]
            stdev = fp0["mleSD",c(3,1,2)]
            pds. = t(ml.out$pds.)
          }
          else if (nulltype==2) {
            pds = fp0["cmest",c(3,1,2)]
            stdev = fp0["cmeSD", c(3,1,2)]
            pds. = t(Cov2.out$pds.)
          }
          else stop("With sw=2, nulltype must equal 0, 1, or 2.")
          colnames(pds.) = names(pds) = c("p0", "delhat", "sighat")
          names(stdev) = c("sdp0", "sddelhat", "sdsighat")
	  return(list(pds=pds, x=x, f=f, pds.=pds., stdev=stdev))
        }
	# find cdf1, the cdf of fdr according to f1 density..................
	p1 <- seq(0.01, 0.99, 0.01)
	cdf1 <- rep(0, 99)
	fd <- fdr
	if(nulltype == 0)
		fd <- fdr0
	for(i in 1:99)
		cdf1[i] <- sum(f1[fd <= p1[i]])
        cdf1 <- cbind(p1, cdf1/cdf1[99])
	mat <- cbind(x, fdr, Fdrl, Fdrr, f, f0, f00, fdr0, yall, lfdrse,
		f1)
	namat <- c("x", "fdr", "Fdrleft", "Fdrright", "f", "f0", "f0theo",
		"fdrtheo", "counts", "lfdrse", "p1f1")
	if(nulltype == 0)
		namat[c(3, 4, 10)] <- c("Fdrltheo", "Fdrrtheo", 
			"lfdrsetheo")
	dimnames(mat) <- list(NULL, namat)
        ############## Locations of triangles ###########
        z.2 = rep(NA, 2)
        m = order(fd)[nx]
        if (fd[nx] < 0.2) {
          z.2[2] = approx(fd[m:nx], x[m:nx], 0.20000000000000001,
                               ties=mean)$y
        }
        if (fd[1] < 0.2) {
          z.2[1] = approx(fd[1:m], x[1:m], 0.20000000000000001,
                              ties=mean)$y
        }
        ################### Plotting ####################
	if(plot > 0) {
		if(plot == 2 | plot == 3)
			oldpar <- par(mfrow = c(1, 2), pty = "m")
                else if (plot ==4) oldpar = par(mfrow = c(1, 3), pty = "m")
		hist(zzz, breaks = breaks, xlab = " ", main = main)
		################### make yt positive ##############
		yt <- pmax(yall * (1 - fd), 0)
		for(k in 1:K)
		  lines(c(x[k], x[k]), c(0, yt[k]), lwd = 2, col = 6)
		if(nulltype == 3)
			title(xlab = paste("delta=", round(xmax, 3),
				"sigleft=", round(sigs[1], 3), 
				" sigright=", round(sigs[2], 3), "p0=",
				round(fp0["cmest", 3], 3)))
		if(nulltype == 1 | nulltype == 2)
			title(xlab = paste("MLE: delta:", round(
				mlests[1], 3), "sigma:", round(mlests[
				2], 3), "p0:", round(mlests[3], 3)),
                                sub = paste("CME: delta:",
                                            round(fp0["cmest",1], 3),
                                            "sigma:", round(fp0["cmest",2], 3),
                                            "p0:", round(fp0["cmest", 3], 3)))
		lines(x, f, lwd = 3, col = 3)
		if(nulltype == 0)
			lines(x, p0theo*f00, lwd = 2, lty = 2, col = 4)
		else
			lines(x, p0*f0, lwd = 2, lty = 2, col = 4)
                ################## Plot triangles ###############
                if (!is.na(z.2[2]))
		   points(z.2[2], -0.5, pch = 24, col="red", bg="yellow")
                if(!is.na(z.2[1]))
		   points(z.2[1], -0.5, pch = 24, col="red", bg="yellow")
		if(nulltype == 1 | nulltype ==2)
			Ef <- Efdr[1]
		else if (nulltype == 0) Ef <- Efdr[4]
		if (plot == 2 | plot == 4) {
			if(nulltype == 0)
				fdd <- fdr0
			else fdd = fdr
			matplot(x, cbind(fdd, Fdrl, Fdrr), type = "l",
				lwd = 3, xlab = " ", ylim = c(0, 
				1.1000000000000001), main = 
				"fdr (solid); Fdr's (dashed)")
			title(xlab = paste("Efdr= ", round(Ef, 3)))
			abline(0, 0, lty = 3, col = 2)
			lines(c(0, 0), c(0, 1), lty = 3, col = 2)
		}
		if (plot == 3 | plot == 4) {
			if(sum(is.na(cdf1[, 2])) == nrow(cdf1))
				warning("cdf1 not available")
			else {
				plot(cdf1[, 1], cdf1[, 2], type = "l",
					lwd = 3, xlab = "fdr level", ylim
					 = c(0, 1), ylab = 
					"f1 proportion < fdr level", main
					 = "f1 cdf of estimated fdr")
				title(sub = paste("Efdr= ", round(Ef,
					3)))
				lines(c(0.20000000000000001, 
					0.20000000000000001), c(0, cdf1[
					20, 2]), col = 4, lty = 2)
				lines(c(0, 0.20000000000000001), rep(
					cdf1[20, 2], 2), col = 4, lty = 2)
				text(0.050000000000000003, cdf1[20, 2],
					round(cdf1[20, 2], 2))
				abline(0, 0, col = 2)
				lines(c(0, 0), c(0, 1), col = 2)
			}
		}
		if(plot > 1)
			par(oldpar)
	}
	if(nulltype == 0) {
		ffdr <- approx(x, fdr0, zz, rule = 2, ties="ordered")$y
	}
	else ffdr <- approx(x, fdr, zz, rule = 2, ties="ordered")$y
	vl = list(fdr = ffdr, fp0 = fp0, Efdr = Efdr, cdf1 = cdf1,  mat = mat,
          z.2 = z.2)
	if(!missing(mult))
		vl$mult = EE
        vl$call = call
	vl
}

