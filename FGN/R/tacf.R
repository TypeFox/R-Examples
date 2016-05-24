mix <- function(x, y) {
	
	n <- 2*length(x)-2
	rev(Re(fft(fft(symtacvf(x)) * fft(symtacvf(y)), inverse = TRUE)/n)[(n/2 - 1):(n - 1)])

}


"tacvfARFIMA" <- 
function(phi = numeric(0), theta = numeric(0), dfrac = numeric(0), phiseas = numeric(0), thetaseas = numeric(0), 
	dfs = numeric(0), H = numeric(0), Hs = numeric(0), alpha = numeric(0), alphas = numeric(0), period = 0, maxlag, useCt = T, sigma2 = 1)
{
	if(length(maxlag) == 0) stop("maxlag must be defined")
	#if(!IdentInvertQ(phi = phi, theta = theta, dfrac = dfrac, phiseas = phiseas, thetaseas = thetaseas, dfs = dfs, H = H, Hs = Hs, alpha = alpha, alphas = alphas, period = period, ident = F)) {
	#	return(NULL) 
	#}

	if(length(period)==0||is.na(period)||period == 0) {
		if(length(phiseas) > 0 && phiseas != 0) stop("Period = 0, but seasonal phi != 0")
		if(length(thetaseas) > 0 && thetaseas != 0) stop("Period = 0, but seasonal theta != 0")
		if((length(dfs) > 0 && dfs != 0)||(length(Hs)>0&&Hs!=0)||(length(alphas)>0&&alphas!=0)) stop("Period = 0, but fractional seasonal parameter != 0")
		return(tacvfFARMA(phi = phi, theta = theta, dfrac = dfrac, H=H, alpha = alpha, maxlag = maxlag, useCt = useCt, sigma2 = sigma2))
	}
	if((period != round(period))||(period < 2)) stop("Period must be an integer >= 2")
	lagTrunc <- 2*max(128, nextn(maxlag, factors = 2))
	iseas <- nextn(period, factors = 2)
	model <- tacvfFARMA(phi = phi, theta = theta, dfrac = dfrac, H=H, alpha = alpha, maxlag = (lagTrunc*iseas), nolagtrunc = T, useCt = useCt, sigma2 = 1)
	modelseas <- tacvfFARMA(phi = phiseas, theta = thetaseas, dfrac = dfs, H = Hs, alpha = alphas, maxlag = (lagTrunc*2), nolagtrunc = T, useCt = useCt, sigma2 = 1)

	modelseas <- shift(modelseas, period, useCt = T) 
	if(length(modelseas) < length(model)) stop("oops")
	modelseas <- modelseas[1:(lagTrunc*iseas+1)]

	z <- sigma2*mix(model, modelseas)
	return(z[1:(maxlag+1)])
}


"tacvfFARMA"   <-  
function(phi = numeric(0), theta = numeric(0), dfrac = numeric(0), H = numeric(0), alpha = numeric(0), maxlag, nolagtrunc = F, useCt = T, sigma2 = 1)
{
	if(length(H)+length(dfrac)+length(alpha)>1) stop("more than one GHD process specified:  stopping")
	if((length(H)==0)&&(length(dfrac) == 0) && (length(alpha) == 0) && (length(phi) == 0) && (length(theta) == 0)) {
		return(c(sigma2, rep(0, maxlag)))
	}
	
	if((length(dfrac) == 0)&&(length(H) == 0)&&(length(alpha) == 0)) {
		return(tacvfARMA(phi = phi, theta = theta, maxlag = maxlag, useCt = useCt, sigma2 = sigma2))
	}
	
	onlyHD <- (length(phi) == 0 && length(theta) == 0)
	
	if(nolagtrunc || onlyHD) lagTrunc <- maxlag
	
	else {
		lagTrunc <- 2*max(128, nextn(maxlag, factors = 2))
	}
	
	
	if(length(dfrac)> 0) x <- tacvfFDWN(dfrac = dfrac, maxlag = lagTrunc, useCt = useCt)
	else if(length(H)>0) x <- tacvfFGN(H=H, maxlag = lagTrunc, useCt = useCt)
	else x <- tacvfHD(alpha = alpha, maxlag = lagTrunc, useCt = useCt)
	
	if(onlyHD) return(sigma2*x[1:(maxlag+1)])
	
	y  <-  tacvfARMA(phi = phi, theta = theta, maxlag = lagTrunc, useCt = useCt, sigma2 = 1)
	
	z <- sigma2*mix(x, y)
	return(z[1:(maxlag+1)])
}

"symtacvf"   <-  
function(x)
{
	c(rev(x[-1])[-1], x)
}

"tacvfARMA"   <-  
function(phi = numeric(0), theta = numeric(0), maxlag = 20, useCt = T, sigma2 = 1)
{

	if(useCt) {
		dd <- function(x, y, ml) .C("tacvfARMA_C", as.double(x), as.integer(length(x)), as.double(y), as.integer(length(y)), as.integer(ml), res = double(ml+1))$res
		return(sigma2*dd(phi, theta, maxlag))
	}
	p   <-   length(phi)
	q   <-   length(theta)
	maxlagp1   <-   maxlag + 1

	if(max(p, q) == 0) {
			return(c(sigma2, numeric(maxlag)))
	}
	r   <-   max(p, q) + 1
	b   <-   numeric(r)
	C   <-   numeric(q + 1)

	C[1]   <-   1
	theta2   <-   c(-1, theta)
	phi2   <-   numeric(3 * r)
	phi2[r]   <-   -1
	if(p > 0) {
		phi2[r + 1:p]   <-   phi
	}
	if(q > 0) {
		for(k in 1:q) {
			C[k + 1]   <-    - theta[k]
			if(p > 0) {
				for(i in 1:min(p, k)) {
				  C[k + 1]   <-   C[k + 1] + phi[i] * C[k + 1 - i]
				}
			}
		}
	}

	for(k in 0:q) {
		for(i in k:q) {
			b[k + 1]   <-   b[k + 1] - theta2[i + 1] * C[i - k + 1]
		}
	}

	if(p == 0) {
		g   <-   c(b, numeric(maxlagp1))[1:maxlagp1]
		
		return(g)
	}
	else if(p > 0) {
		a   <-   matrix(numeric(r^2), ncol = r)
		for(i in 1:r) {
			for(j in 1:r) {
				if(j == 1) {
				  a[i, j]   <-   phi2[r + i - 1]
				}
				else if(j != 1) {
				  a[i, j]   <-   phi2[r + i - j] + phi2[r + i + j - 
				    2]
				}
			}
		}


		g   <-   solve(a,  - b)
		
		if(length(g) <= maxlag) {
			g   <-   c(g, numeric(maxlagp1 - r)) 

			for(i in (r + 1):maxlagp1) {
				g[i]   <-   phi %*% g[i - 1:p]
			}
		
			return(sigma2*g[1:maxlagp1])
		}
		else if(length(g) >= maxlagp1) {
			
			return(sigma2*g[1:maxlagp1])
		}
	}
}

"tacvfFDWN"   <-  
function(dfrac, maxlag, useCt = T)
{

	#if(!InvertibleD(dfrac)) {
	#	warning("Model is non-causal or non-invertible\n")
	#		return(NULL)
	#}
	if(useCt) {
		ta <- function(ds, ma) .C("tacvfFDWN_C", as.double(ds), as.integer(ma), x = double(ma+1))$x
		return(ta(dfrac, maxlag))
	}
	x   <-   numeric(maxlag + 1)
	x[1]   <-   gamma(1 - 2 * dfrac)/gamma(1 - dfrac)^2
	for(i in 1:maxlag) {
		x[i + 1]   <-   ((i - 1 + dfrac)/(i - dfrac)) * x[i]
	}
	return(x)
}

tacvfFGN <-
function(H, maxlag, useCt = TRUE){
	#if(!InvertibleH(H)) return(NULL)
	if(useCt) {
		tg <- function(H, ma) .C("tacfFGN_C", as.double(H), as.integer(ma), x = double(ma+1))$x
		return(tg(H, maxlag))
	}
	h2<-2*H 
    r <- sapply(0:maxlag, function(k) 0.5*(abs(k+1)^h2-2*abs(k)^h2+abs(k-1)^h2))
	return(r)
}
 
## will lgamma make it faster??
Zeta <- function(s, n = 20) {

	d <- rep(0, n+1)
	
	d[1] <- n * (factorial(n - 1))/(factorial(n))
	
	for(i in 1:(n)) {
		d[i+1] <- d[i] + n * factorial(n + i -1)*4^i/(factorial(n-i)*factorial(2*i))
	}

	zeta <- 0
	
	for(i in 0:(n-1)) {
		zeta <- zeta + (-1)^i*(d[i+1] - d[n+1])/(i+1)^s
	}
	
	zeta <- zeta * (-1)/(d[n+1]*(1 - 2^(1-s)))
	zeta
}
##Reference paper. (zetaalgorithm.)

tacvfHD <- 
function(alpha, maxlag, useCt = TRUE) {
	#if(!InvertibleAlpha(alpha)) return(NULL)
	if(useCt) {
		th <- function(alpha, ma) .C("tacfHD_C", as.double(alpha), as.integer(ma), x = double(ma+1))$x
		return(th(alpha, maxlag))
	}
	val <-  (-2*Zeta(alpha))^(-1)
	r <- c(1, val*sapply(1:maxlag, function(k) k^(-alpha)))
	return(r)
}


psiwts <- function(phi = 0, theta = 0, phiseas = 0, thetaseas = 0, dfrac = 0, dfs = 0, dint = 0, dseas = 0, period = 0, len = 128, n = len*2, div = 1)
{
	multer = 1
	if((length(theta) > 0) && any(theta != 0)) multer <- mult(multer, c(1, -theta))
	if((length(phi) > 0) && any(phi != 0)) multer <- mult(multer, psiwtsAR(phi, n))
	if(length(dint) > 0 && dint > 0 && length(dfrac) > 0 && dfrac != 0) 
		multer <- mult(multer, expand(d = -(dint+dfrac), n = n))
	else if(length(dint) > 0 && dint > 0)
		multer <- mult(multer, expand(d = -dint, n = n))
	else if(length(dfrac) > 0 && dfrac != 0)
		multer <- mult(multer, expand(d = -dfrac, n = n))

	if(period > 0) {
		if(length(phiseas) > 0 && any(phiseas != 0)) multer <- mult(multer, psiwtsAR(shift(c(1, phiseas), period)[-1], n*period))
		if(length(thetaseas) > 0 && any(thetaseas != 0)) multer <- mult(multer, shift(c(1, -thetaseas), period))
		if(length(dseas) > 0 && dseas > 0 && length(dfs) > 0 && dfs != 0) 
			multer <- mult(multer, expandseas(d = -(dseas+dfs), seas = period, n = n))
		else if(length(dseas) > 0 && dseas > 0)
			multer <- mult(multer, expandseas(d = -dseas, seas = period, n = ceiling(n/div)))
		else if(length(dfs) > 0 && dfs != 0)
			multer <- mult(multer, expandseas(d = -dfs, seas = period, n = n))
	}
	if(length(multer) < len) multer <- c(multer, rep(0, len))
	return(multer[1:len])
}

"psiwtsAR" <-
function(phi, maxlag)
{
	p <- length(phi)
	if(p == 0) return(0)
	x <- numeric(maxlag + 1)
	x <- 1
	for(i in 1:p) {
		x[i + 1] <- crossprod(phi[1:i], (rev(x))[1:i])
	}
	if(maxlag > p) {
		for(i in (p + 1):maxlag) {
			x[i + 1] <- crossprod(phi, (rev(x))[1:p])
		}
	}
	return(x)
}

"shift" <- function(pars, by.am, useCt = TRUE) {
	
	if(useCt && by.am > 0) {
		ss <- function(vv, byam)  .C("shift_C", as.double(vv), as.integer(byam), as.integer(length(vv)), y = double((length(vv)-1)*byam+1))$y
		return(ss(pars, by.am))
	}
	n <- length(pars)
	
	if(all(pars[-1] == 0))
		return(pars)
	
	if(by.am < 0) {
		by.am <- abs(by.am)
		m <- floor((n-1)/by.am)
		
		seq.pol <- rep(0, m+1)
		
		seq.pol[1] <- pars[1]
		
		for(i in 2:(m+1)) seq.pol[i] <- pars[(i-1)*by.am+1]
		
		return(seq.pol)
		
	}
	
	seq.pol <- rep(0, (n-1)*by.am+1)
	
	seq.pol[1] <- pars[1]
	
	for(i in 2:n) seq.pol[(i-1)*by.am+1] <- pars[i]
	
	return(seq.pol)
}

"mult" <- function(a, b) {
	if(all(a == 0)||all(b == 0)) return(0)
	ret <- convolve(a, rev(b), type = "open")
	return(ret)
}

expand = function(d, Bterm = -1, n = 10) {
	
	ret = as.vector(sapply(0:n, function(x) choose(d, x)*Bterm^x))
	
	return(ret)
	
}

expandseas = function(d, Bterm = -1, seas = 12, n = 10) {
	
	return(shift(expand(d, Bterm, n), seas))
	
}
