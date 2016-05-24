# This is ITSM-R.
#
# Notes:
#
# 1. The acronym ITSF appears in many of the comments below. It stands for
#    "Introduction to Time Series and Forecasting" by Brockwell and Davis.
#
# 2. The ITSM-R forecast function is designed to replicate the results of the
#    Windows program B&D ITSM. These results will differ from the standard R
#    function predict.arima().

selftest = function() {
	.test.smooth.ma()
	.test.smooth.exp()
	.test.smooth.fft()
	.test.acvf()
	.test.aacvf()
	.test.ma.inf()
	.test.ar.inf()
	.test.ia()
	.test.forecast.lake()
	.test.forecast.dowj()
	.test.forecast.wine()
	.test.forecast.deaths()
	.test.forecast.sunspots()
	.test.arar.deaths()
	.test.arar.trings()
	.test.periodogram()
}

wine = c(
464,675,703,887,1139,1077,1318,1260,1120,963,996,960,530,883,894,1045,
1199,1287,1565,1577,1076,918,1008,1063,544,635,804,980,1018,1064,1404,
1286,1104,999,996,1015,615,722,832,977,1270,1437,1520,1708,1151,934,
1159,1209,699,830,996,1124,1458,1270,1753,2258,1208,1241,1265,1828,809,
997,1164,1205,1538,1513,1378,2083,1357,1536,1526,1376,779,1005,1193,
1522,1539,1546,2116,2326,1596,1356,1553,1613,814,1150,1225,1691,1759,
1754,2100,2062,2012,1897,1964,2186,966,1549,1538,1612,2078,2137,2907,
2249,1883,1739,1828,1868,1138,1430,1809,1763,2200,2067,2503,2141,2103,
1972,2181,2344,970,1199,1718,1683,2025,2051,2439,2353,2230,1852,2147,
2286,1007,1665,1642,1525,1838,1892,2920,2572,2617,2047)

deaths = c(
9007,8106,8928,9137,10017,10826,11317,10744,9713,9938,9161,8927,7750,
6981,8038,8422,8714,9512,10120,9823,8743,9129,8710,8680,8162,7306,8124,
7870,9387,9556,10093,9620,8285,8433,8160,8034,7717,7461,7776,7925,8634,
8945,10078,9179,8037,8488,7874,8647,7792,6957,7726,8106,8890,9299,10625,
9302,8314,8850,8265,8796,7836,6892,7791,8129,9115,9434,10484,9827,9110,
9070,8633,9240)

strikes = c(
4737,5117,5091,3468,4320,3825,3673,3694,3708,3333,3367,3614,3362,3655,
3963,4405,4595,5045,5700,5716,5138,5010,5353,6074,5031,5648,5506,4230,
4827,3885)

lake = c(
10.38,11.86,10.97,10.8,9.79,10.39,10.42,10.82,11.4,11.32,11.44,11.68,
11.17,10.53,10.01,9.91,9.14,9.16,9.55,9.67,8.44,8.24,9.1,9.09,9.35,8.82,
9.32,9.01,9,9.8,9.83,9.72,9.89,10.01,9.37,8.69,8.19,8.67,9.55,8.92,8.09,
9.37,10.13,10.14,9.51,9.24,8.66,8.86,8.05,7.79,6.75,6.75,7.82,8.64,
10.58,9.48,7.38,6.9,6.94,6.24,6.84,6.85,6.9,7.79,8.18,7.51,7.23,8.42,
9.61,9.05,9.26,9.22,9.38,9.1,7.95,8.12,9.75,10.85,10.41,9.96,9.61,8.76,
8.18,7.21,7.13,9.1,8.25,7.91,6.89,5.96,6.8,7.68,8.38,8.52,9.74,9.31,
9.89,9.96)

Sunspots = c(
101,82,66,35,31,7,20,92,154,125,85,68,38,23,10,24,83,132,131,118,90,67,
60,47,41,21,16,6,4,7,14,34,45,43,48,42,28,10,8,2,0,1,5,12,14,35,46,41,
30,24,16,7,4,2,8,17,36,50,62,67,71,48,28,8,13,57,122,138,103,86,63,37,
24,11,15,40,62,98,124,96,66,64,54,39,21,7,4,23,55,94,96,77,59,44,47,30,
16,7,37,74)

dowj = c(
110.94,110.69,110.43,110.56,110.75,110.84,110.46,110.56,110.46,110.05,
109.6,109.31,109.31,109.25,109.02,108.54,108.77,109.02,109.44,109.38,
109.53,109.89,110.56,110.56,110.72,111.23,111.48,111.58,111.9,112.19,
112.06,111.96,111.68,111.36,111.42,112,112.22,112.7,113.15,114.36,
114.65,115.06,115.86,116.4,116.44,116.88,118.07,118.51,119.28,119.79,
119.7,119.28,119.66,120.14,120.97,121.13,121.55,121.96,122.26,123.79,
124.11,124.14,123.37,123.02,122.86,123.02,123.11,123.05,123.05,122.83,
123.18,122.67,122.73,122.86,122.67,122.09,122,121.23)

airpass = c(
112,118,132,129,121,135,148,148,136,119,104,118,115,126,141,135,125,149,
170,170,158,133,114,140,145,150,178,163,172,178,199,199,184,162,146,166,
171,180,193,181,183,218,230,242,209,191,172,194,196,196,236,235,229,243,
264,272,237,211,180,201,204,188,235,227,234,264,302,293,259,229,203,229,
242,233,267,269,270,315,364,347,312,274,237,278,284,277,317,313,318,374,
413,405,355,306,271,306,315,301,356,348,355,422,465,467,404,347,305,336,
340,318,362,348,363,435,491,505,404,359,310,337,360,342,406,396,420,472,
548,559,463,407,362,405,417,391,419,461,472,535,622,606,508,461,390,432)

.trings = c(
1.046,1.038,0.925,0.755,1.156,1.163,0.933,1.127,1.053,1.191,0.963,0.961,
0.729,0.898,0.925,0.802,0.776,0.865,0.87,0.906,0.818,0.76,0.709,0.642,
0.797,0.585,0.782,0.587,0.77,0.512,0.702,0.795,1.088,1.074,1.06,1,1.091,
1.209,1.077,0.481,0.974,1.219,1.213,1.05,1.05,0.92,1.32,1.033,1.309,
0.979,1.15,0.908,1.113,0.962,1.283,1.094,1.08,0.956,1.104,0.984,1.028,
1.159,1.223,0.933,1.167,1.114,0.923,1.03,1.126,0.87,0.824,1.277,1.182,
1.365,1.269,1.23,1.446,1.602,1.328,1.294,1.208,1.203,1.51,1.191,1.069,
0.927,0.828,0.98,1.187,1.307,1.003,1.128,1.198,1.148,1.19,1.024,1.272,
1.283,1.225,1.269,1.14,0.963,0.713,0.817,1.401,1.146,0.826,0.812,0.866,
0.997,1.091,1.157,1.077,0.951,1.061,0.899,0.817,0.922,0.579,0.863,0.697,
0.814,0.835,0.938,1.201,1.261,0.71,0.6,0.576,0.606,0.579,0.765,0.823,
0.856,1.041,0.888,0.674,0.73,0.81,0.73,0.804,0.772,0.741,0.955,0.694,
0.829,0.83,0.89,0.884,0.923,0.865,1.026,1.043,0.537,0.693,0.904,0.902,
0.779,0.759,0.707,0.804,0.763,0.517,0.704,0.584,0.372,0.607,0.81,0.733,
0.618,0.65,0.603,0.702,0.9,0.883,0.727,0.812,0.673,0.724,0.958,0.954,
0.845,0.638,0.987,0.987,1.098,1.287,1.179,1.193,1.246,1.084,1.141,1.245,
1.174,1.054,0.924,1.118,1.149,1.108,1.305,1.314,1.279,1.106,1.145,0.977,
1.427,1.13,1.152,1.247,1.344,1.22,1.198,1.253,0.919,1.057,1.231,1.286,
1.354,1.277,1.553,1.123,1.327,0.525,0.966,0.952,0.778,0.765,1.328,0.753,
1.243,0.943,1.257,1.122,1.265,1.367,1.247,1.808,2.02,1.383,1.16,1.017,
0.976,0.771,0.662,0.946,1.148,0.948,1.026,1.248,0.949)

# Compute the autocovariance vector for an ARMA model
#
# Solve first p+1 equations to obtain gamma(0), ..., gamma(p).
#
# Example: p=2
#
# k=0: gamma(0) - phi1 * gamma(-1) - phi2 * gamma(-2) = rhs(0)
#
# k=1: gamma(1) - phi1 * gamma(0)  - phi2 * gamma(-1) = rhs(1)
#
# k=2: gamma(2) - phi1 * gamma(1)  - phi2 * gamma(0)  = rhs(2)
#
# Hence
#
#	gamma(2)	gamma(1)	gamma(0)	gamma(-1)	gamma(-2)
#
# k=0:	0		0		1		-phi1		-phi2
#
# k=1:	0		1		-phi1		-phi2		0
#
# k=2:	1		-phi1		-phi2		0		0
#
# Since gamma is symmetric, fold the matrix in half by adding columns:
#
#					gamma(0)	gamma(-1)	gamma(-2)
#
# k=0:					1		-phi1		-phi2
#
# k=1:					-phi1		1-phi2		0
#
# k=2:					-phi2		-phi1		1
#
# Arguments
#
#	a	List of $phi, $theta, $sigma2
#
#	h	Maximum lag
#
# Returns a vector of length h+1 to accommodate lag 0 at index 1.

aacvf = function(a,h) {

	phi = a$phi
	theta = a$theta
	sigma2 = a$sigma2

	p = length(phi)
	q = length(theta)

	psi = ma.inf(a,q)
	theta = c(1,theta)

	# r is the right hand side of equation (3.2.5)

	f = function(k) sum(theta[k:(q+1)]*psi[1:(q-k+2)])
	r = numeric(max(p+1,q+1,h+1))
	r[1:(q+1)] = sigma2*sapply(1:(q+1),f)

	# Solve for gamma in A * gamma = r

	f = function(k) c(numeric(p-k),1,-phi,numeric(k))
	A = sapply(0:p,f)
	A = matrix(A,ncol=2*p+1,byrow=TRUE)
	A = cbind(A[,p+1],A[,(p+2):(2*p+1)]+A[,p:1]) #Fold A

	gamma = numeric(max(p+1,h+1))
	gamma[1:(p+1)] = solve.qr(qr(A),r[1:(p+1)])

	# Calculate remaining lags recursively

	if (h > p)
		for (k in (p+1):h)
			 gamma[k+1] = r[k+1] + sum(phi*gamma[k-(1:p)+1])
	else if (h < p)
		gamma = gamma[1:(h+1)]

	return(gamma)
}

# Test aacvf() using ARMAacf()

.test.aacvf = function() {
	cat("aacvf ")
	for (p in 0:6) {
		for (q in 0:6) {
			if (p == 0)
				phi = 0
			else
				phi = 0.1 + (1:p)/100
			if (q == 0)
				theta = 0
			else
				theta = 0.2 + (1:q)/100
			a = list(phi=phi,theta=theta,sigma2=1)
			b = aacvf(a,40)
			b = b/b[1]
			c = ARMAacf(ar=phi,ma=theta,lag.max=40)
			if (any(abs(b-c) > 1e-9)) {
				cat("fail",p,q,"\n")
				return()
			}
		}
	}
	cat("ok\n")
}

# Returns autocovariance of data
#
# Arguments
#
#	x	Data
#
#	h	Maximum lag
#
# Returns a vector of length h+1 to accommodate lag 0 at index 1.

acvf = function(x,h=40) {
	n = length(x)
	if (h >= n)
		stop("lag > data")
	xbar = mean(x)
	f = function(i) sum((x[1:(n-i)]-xbar) * (x[(i+1):n]-xbar)) / n
	return(sapply(0:h,f))
}

# Test acvf() using acf()

.test.acvf = function() {
	cat("acvf ")
	a = acvf(Sunspots,40)
	b = acf(Sunspots,lag.max=40,type="covariance",plot=FALSE)$acf
	if (any(abs(a-b)>1e-9)) cat("fail\n") else cat("ok\n");
}

# Plot one or two lines in color similar to ITSM
#
# Arguments
#
#	y1	Blue line
#
#	y2	Red line

plotc = function(y1,y2=NULL) {
	if (is.null(y2)) {
		if (is.ts(y1))
			x1 = time(y1)
		else
			x1 = 1:length(y1)
		plot.default(x1,y1,xlab="",ylab="",type="o",col="blue")
		return(invisible(NULL))
	}
	n1 = length(y1)
	n2 = length(y2)
	if (is.ts(y1)) {
		x1=time(y1)
		x2=(0:(n2-1))/frequency(y1)+start(y1)
	} else {
		x1=seq(1,n1)
		x2=seq(1,n2)
	}
	xmax = max(x1,x2)
	xmin = min(x1,x2)
	ymax = max(y1,y2)
	ymin = min(y1,y2)
	plot.default(x1,y1,xlim=c(xmin,xmax),ylim=c(ymin,ymax),xlab="",ylab="",type="o",col="blue")
	lines(x2,y2,type="l",col="red")
}

# Plot data and forecasted values
#
# Arguments
#
#	x	Observed data
#
#	f	List of $pred, $l, $u
#
# If x is a time series object then time(x) labels are used.

.plot.forecast = function(x,f) {
	n = length(x)
	m = length(f$pred)
	ymin = min(x,f$l)
	ymax = max(x,f$u)
	if (is.ts(x)) {
		u = time(x)
		v = (0:(n+m-1))/frequency(x)+start(x)
	} else {
		u = 1:n
		v = 1:(n+m)
	}
	plot.default(u,x,xlim=c(v[1],v[n+m]),ylim=c(ymin,ymax),xlab="",ylab="",type="o",col="blue")
	lines(v[(n+1):(n+m)],f$pred[1:m],type="o",col="red")
	lines(v[(n+1):(n+m)],f$l[1:m],type="l",col="red",lty=2)
	lines(v[(n+1):(n+m)],f$u[1:m],type="l",col="red",lty=2)
}

# Returns the trend component of data
#
# Arguments
#
#	x	Data
#
#	p	Polynomial order (1 linear, 2 quadratic, etc.)

trend = function(x,p) {
	n = length(x)
	X = NULL
	for (i in (0:p))
		X = cbind(X,(1:n)^i)
	b = solve.qr(qr(X),x)
	xhat = as.vector(X %*% b)
	return(xhat)
}

# Smooth data with a moving average filter
#
# Arguments
#
#	x	Data
#
#	q	Window size is 2q + 1
#
# Return
#
#	Applies a moving average filter to x and returns the result.

smooth.ma = function(x,q) {
	n = length(x)
	x = c(rep(x[1],q),x,rep(x[n],q)) #Extend data per B&D p. 25
	qq = -q:q
	F = function(t) sum(x[t+qq])/(2*q+1)
	m = sapply((q+1):(n+q),F)
	return(m)
}

.test.smooth.ma = function() {
	cat("smooth.ma ")
	d = smooth.ma(strikes,2) - c(
     4883.800000000000000, #Data from B&D ITSM using File->Export
     4630.000000000000000,
     4546.600000000000000,
     4364.200000000001000,
     4075.400000000001000,
     3796.000000000001000,
     3844.000000000000000,
     3646.600000000000000,
     3555.000000000000000,
     3543.200000000000000,
     3476.800000000000000,
     3466.200000000001000,
     3592.200000000000000,
     3799.800000000000000,
     3996.000000000000000,
     4332.600000000000000,
     4741.600000000000000,
     5092.200000000000000,
     5238.800000000000000,
     5321.800000000000000,
     5383.400000000000000,
     5458.200000000001000,
     5321.200000000001000,
     5423.200000000001000,
     5522.400000000001000,
     5297.800000000000000,
     5048.400000000001000,
     4819.200000000000000,
     4466.600000000000000,
     4142.400000000000000)
	if (any(abs(d) > 1e-6))
		stop("fail\n")
	cat("ok\n")
}

# Exponential smoothing
#
# Arguments
#
#	x	Data
#
#	alpha	Smoothing parameter 0-1
#
# Return
#
#	Applies an exponential filter to x and returns the result.

smooth.exp = function(x,alpha) {
	n = length(x)
	m = numeric(n)
	m[1] = x[1]
	for (t in 2:n) #Recursive algorithm, for-loop required
		m[t] = alpha*x[t]+(1-alpha)*m[t-1]
	return(m)
}

.test.smooth.exp = function() {
	cat("smooth.exp ")
	d = smooth.exp(strikes,0.4) - c(
     4737.000000000000000, #Data from B&D ITSM using File->Export
     4889.000000000000000,
     4969.800000000000000,
     4369.080000000000000,
     4349.448000000000000,
     4139.668800000000000,
     3953.001280000000000,
     3849.400768000000000,
     3792.840460800000000,
     3608.904276480000000,
     3512.142565888000000,
     3552.885539532800000,
     3476.531323719680000,
     3547.918794231808000,
     3713.951276539085000,
     3990.370765923451000,
     4232.222459554070000,
     4557.333475732442000,
     5014.400085439465000,
     5295.040051263679000,
     5232.224030758207000,
     5143.334418454924000,
     5227.200651072955000,
     5565.920390643772000,
     5351.952234386264000,
     5470.371340631758000,
     5484.622804379055000,
     4982.773682627433000,
     4920.464209576459000,
     4506.278525745875000)
	if (any(abs(d) > 1e-6))
		stop("fail\n")
	cat("ok\n")
}

# FFT smoothing
#
# Arguments
#
#	x	Data
#
#	f	Fraction of spectrum to pass 0-1
#
# Return
#
#	Applies a low pass filter to x and returns the result.

smooth.fft = function(x,f) {

	# Fourier transform per equation (4.2.5)

	n = length(x)
	w = floor(n/2)
	t = 1:n
	F1 = function(k) sum(x*complex(n,cos(2*pi*k/n*t),-sin(2*pi*k/n*t)))
	a = sapply(-w:w,F1)/sqrt(n)

	# Low-pass filter

	q = floor(f*w)
	a = a * c(rep(0,w-q),rep(1,2*q+1),rep(0,w-q))

	# Inverse Fourier transform per (4.2.6) and divide by sqrt(n)

	k = -w:w
	F2 = function(t) sum(a*complex(2*w+1,cos(2*pi*k/n*t),sin(2*pi*k/n*t)))
	y = sapply(1:n,F2)/sqrt(n)
	y = Re(y)

	return(y)
}

.test.smooth.fft = function() {
	cat("smooth.fft ")
	d = smooth.fft(strikes,0.4) - c(
     4721.322634626735000,
     4927.422792660943000,
     4739.796871252751000,
     4271.865921112883000,
     3866.507707290081000,
     3737.853393151683000,
     3788.164077374752000,
     3775.225836209553000,
     3604.390085643926000,
     3406.621691140404000,
     3349.677238673979000,
     3444.122606820845000,
     3574.924431536113000,
     3686.480095117960000,
     3868.277893566235000,
     4235.239286500041000,
     4757.852560993248000,
     5247.195999898472000,
     5508.024641920165000,
     5495.294566766186000,
     5335.715987959870000,
     5224.027594619088000,
     5288.353152446581000,
     5500.509090438299000,
     5675.730881411944000,
     5589.910042978935000,
     5172.985828460593000,
     4624.078839951079000,
     4293.276006961120000,
     4379.152242806328000)
	if (any(abs(d) > 1e-6))
		stop("fail\n")
	cat("ok\n")
}

# Returns the seasonal component of data
#
# Arguments
#
#	x	Data
#
#	d	Number of observations per season

season = function(x,d) {
	n = length(x)
	q = floor(d/2)
	if (d == 2*q) {
		F1 = function(t) x[t-q]/2 + sum(x[(t-q+1):(t+q-1)]) + x[t+q]/2
		m = sapply((q+1):(n-q),F1) / d
		m = c(rep(0,q),m,rep(0,q))
	} else
		m = smooth.ma(x,q)
	dx = x - m
	F2 = function(k) mean(dx[seq(k+q,n-q,d)])
	w = sapply(1:d,F2)
	w = w - mean(w)
	s = rep(w,len=n+q)[(q+1):(q+n)]
	return(s)
}

# Spectral filter

smooth.rank = function(x,k) {
	n = length(x)
	q = floor(n/2)
	w = fft(x)/n
	o = order(Mod(w[2:(q+1)]))
	o = o[1:(q-k)]
	w[1+o] = 0
	w[n+1-rev(o)] = 0
	Re(fft(w,inverse=TRUE))
}

# Returns the harmonic component of data
#
# Arguments
#
#	x	Data
#
#	d	Vector of cycle periods

hr = function(x,d) {
	n = length(x)
	X = rep(1,n)
	for (i in d)
		X = cbind(X,cos(2*pi*(0:(n-1))/i),sin(2*pi*(0:(n-1))/i))
	b = solve.qr(qr(X[(1:n),]),x)
	xhat = as.vector(X %*% b)
	return(xhat)
}

# MA(infinity)
#
# Arguments
#
#	a	ARMA model
#
#		phi	AR coefficients phi_1, ..., phi_p
#
#		theta	MA coefficients theta_1, ..., theta_q
#
#	n	Required order of psi
#
# Return value
#
#	Coefficient vector of length n+1 to accommodate psi_0 at index 1, where
#
#	psi(z) = theta(z)/phi(z) = psi_0 + psi_1 * z + psi_2 * z^2 + ...
#
# From ITSF p. 85:
#
#                  p
# psi  = theta  + sum phi  psi
#    j        j   k=1    k    j-k
#
# Convert sum to row times column vector:
#
#                                                   T
# psi  = theta  + (phi  ... phi )(psi    ... psi   )
#    j        j       1        p     j-1        j-p
#
# For psi index, convert to sequence plus offset p+1:
#
# ((j-1):(j-p)) + (p+1) = (j+p):(j+1) = (p:1) + j

ma.inf = function(a,n=50) {
	if (n == 0)
		return(1)
	theta = c(a$theta,numeric(n))
	phi = a$phi
	p = length(phi)
	psi = c(numeric(p),1,numeric(n))
	for (j in 1:n)
		psi[j+p+1] = theta[j] + sum(phi * psi[(p:1)+j])
	return(psi[(0:n)+p+1])
}

# Check that ma.inf = theta(z) / phi(z)

.test.ma.inf = function() {
	cat("ma.inf ")
	theta = 0.15
	phi = c(0.1,0.2)
	model = list(phi=phi,theta=theta)
	psi = ma.inf(model,n=30)
	psi.z = function(z) sum(psi * z^(0:30))
	q = length(theta)
	p = length(phi)
	theta.z = function(z) 1 + sum(theta * z^(1:q))
	phi.z = function(z) 1 - sum(phi * z^(1:p))
	z = seq(-1,1,0.1)
	a = sapply(z,psi.z)
	b = sapply(z,theta.z)
	c = sapply(z,phi.z)
	if (any(abs(a-b/c) > 1e-9)) cat("fail\n") else cat("ok\n")
}

# AR(infinity)
#
# Arguments
#
#	a	ARMA model
#
#		phi	AR coefficients phi_1, ..., phi_p
#
#		theta	MA coefficients theta_1, ..., theta_q
#
#	n	Required order of pi
#
# Return value
#
#	Coefficient vector of length n+1 to accommodate pi_0 at index 1, where
#
#	pi(z) = phi(z)/theta(z) = pi_0 + pi_1 * z + pi_2 * z^2 + ...
#
# From ITSF p. 86:
#
#                q
# pi  = -phi  - sum theta  pi
#   j       j   k=1      k   j-k
#
# Convert sum to row times column vector:
#
#                                                   T
# pi  = -phi  - (theta  ... theta )(pi    ... pi   )
#   j       j         1          q    j-1       j-q
#
# For pi index, convert to sequence plus offset q+1:
#
# ((j-1):(j-q)) + (q+1) = (j+q):(j+1) = (q:1) + j

ar.inf = function(a,n=50) {
	if (n == 0)
		return(1)
	phi = c(a$phi,numeric(n))
	theta = a$theta
	q = length(theta)
	pie = c(numeric(q),1,numeric(n))
	for (j in 1:n)
		pie[j+q+1] = -phi[j] - sum(theta * pie[(q:1)+j])
	return(pie[(0:n)+q+1])
}

# Check that ar.inf = phi(z) / theta(z)

.test.ar.inf = function() {
	cat("ar.inf ")
	phi = c(0.12,0.13,0.14)
	theta = c(0.21,0.22)
	model = list("phi"=phi,"theta"=theta)
	pie = ar.inf(model,n=30)
	pie.z = function(z) sum(pie * z^(0:30))
	p = length(phi)
	q = length(theta)
	phi.z = function(z) 1 - sum(phi * z^(1:p))
	theta.z = function(z) 1 + sum(theta * z^(1:q))
	z = seq(-1,1,0.1)
	a = sapply(z,pie.z)
	b = sapply(z,phi.z)
	c = sapply(z,theta.z)
	if (any(abs(a-b/c) > 1e-9)) cat("fail\n") else cat("ok\n")
}

# Plots model and/or data ACF, PACF
#
# Arguments
#
#	u,v	Model and/or data in either order
#
#	h	Maximum lag

plota = function(u,v=NULL,h=40) {

	# Plot data ACF and PACF side by side like ITSM

	plot.data.acf = function(x) {
		n = length(x)
		h = min(h,n-1)
		op = par(mfrow=c(1,2))
		acf = acf(x,lag.max=h,type="correlation",plot=FALSE)$acf
		plot(0:h,acf,ylim=c(-1,1),type="h",xlab="Lag",ylab="",main="ACF",col="blue")
		b = 1.96/sqrt(n)
		abline(h=b,lty=2)
		abline(h=-b,lty=2)
		abline(h=0)
		legend("topright","Data",fill="blue")
		pacf = acf(x,lag.max=h,type="partial",plot=FALSE)$acf
		pacf = c(1,pacf)
		plot(0:h,pacf,ylim=c(-1,1),type="h",xlab="Lag",ylab="",main="PACF",col="blue")
		abline(h=b,lty=2)
		abline(h=-b,lty=2)
		abline(h=0)
		legend("topright","Data",fill="blue")
		par(op)
	}

	# Plot model ACF and PACF side by side like ITSM

	plot.model.acf = function(a) {
		op = par(mfrow=c(1,2))
		acf = ARMAacf(ar=a$phi,ma=a$theta,lag.max=h)
		pacf = ARMAacf(ar=a$phi,ma=a$theta,lag.max=h,pacf=TRUE)
		pacf = c(1,pacf)
		plot(0:h,acf,ylim=c(-1,1),type="h",xlab="Lag",ylab="",main="ACF",col="red")
		abline(h=0)
		legend("topright","Model",fill="red")
		plot(0:h,pacf,ylim=c(-1,1),type="h",xlab="Lag",ylab="",main="PACF",col="red")
		abline(h=0)
		legend("topright","Model",fill="red")
		par(op)
	}

	# Plot both data and model ACF and PACF side by side like ITSM

	plot.combo.acf = function(x,a) {
		n = length(x)
		h = min(h,n-1)
		op = par(mfrow=c(1,2))
		acf = acf(x,lag.max=h,type="correlation",plot=FALSE)$acf
		plot(0:h,acf,ylim=c(-1,1),type="h",xlab="Lag",ylab="",main="ACF",col="blue")
		acf = ARMAacf(ar=a$phi,ma=a$theta,lag.max=h)
		for (i in 0:h) lines(c(i+0.2,i+0.2),c(0,acf[i+1]),col="red")
		b = 1.96/sqrt(n)
		abline(h=b,lty=2)
		abline(h=-b,lty=2)
		abline(h=0)
		legend("topright",c("Data","Model"),fill=c("blue","red"))
		pacf = acf(x,lag.max=h,type="partial",plot=FALSE)$acf
		pacf = c(1,pacf)
		plot(0:h,pacf,ylim=c(-1,1),type="h",xlab="Lag",ylab="",main="PACF",col="blue")
		pacf = ARMAacf(ar=a$phi,ma=a$theta,lag.max=h,pacf=TRUE)
		pacf = c(1,pacf)
		for (i in 0:h) lines(c(i+0.2,i+0.2),c(0,pacf[i+1]),col="red")
		abline(h=b,lty=2)
		abline(h=-b,lty=2)
		abline(h=0)
		legend("topright",c("Data","Model"),fill=c("blue","red"))
		par(op)
	}

	if (is.null(v))
		if (is.list(u))
			plot.model.acf(u)
		else
			plot.data.acf(u)
	else
		if (is.list(u))
			plot.combo.acf(v,u)
		else
			plot.combo.acf(u,v)
}

# Plot data or ARMA model spectrum
#
# Argument
#
#	u	Data or ARMA model (list of $phi, $theta)

plots = function(u) {

	plot.data.spectrum = function(x) {
		n = length(x)
		q = floor(n/2)
		A = Mod(fft(x)/n)
		plot(0:q,A[1:(q+1)],type="h",col="blue",xlab="",ylab="",main="Data Spectrum")
	}

	plot.model.spectrum = function(a) {
		h = 1000
		gamma = aacvf(a,h)
		gamma0 = gamma[1]
		gamma = gamma[2:(h+1)]
		k = 1:h
		f = function(lambda) (gamma0 + 2*sum(cos(k*lambda)*gamma))/2/pi
		lambda = pi * seq(0,1,0.001)
		y = sapply(lambda,f)
		plot(lambda,y,type="l",xlab="",ylab="",col="red",main="Model Spectrum")
	}

	if (is.list(u))
		plot.model.spectrum(u)
	else
		plot.data.spectrum(u)
}

# Specify an ARMA model
#
# Arguments
#
#	ar	AR coefficients [1:p] = phi_1, ..., phi_p
#
#	ma	MA coefficients [1:q] = theta_1, ..., theta_q
#
#	sigma2	White noise variance

specify = function(ar=0,ma=0,sigma2=1) {
	a = list(phi=ar,theta=ma,sigma2=sigma2)
	return(a)
}

# Returns synthetic observations for an ARMA process
#
# This is a wrapper for arima.sim()
#
# Arguments
#
#	a	ARMA model (list of $phi, $theta, $sigma2)
#
#	n	Number of observations required

sim = function(a,n=100) {
	if (all(a$theta == 0))
		y = rnorm(1000+n,mean=0,sd=sqrt(a$sigma2))
	else {
		y = rnorm(1100+n,mean=0,sd=sqrt(a$sigma2))
		y = filter(y,c(1,a$theta),sides=1)
		y = y[-(1:100)]
	}
	if (any(a$phi != 0))
		y = filter(y,a$phi,method="recursive")
	return(y[-(1:1000)]) ## remove first 1000 values
}

# Estimate MA coeffcients using the innovations algorithm
#
# Arguments
#
#	x	Data
#
#	q	MA order
#
#	m	Recursion level
#
# Return value
#
#	$phi	0
#
#	$theta	MA coefficients
#
#	$sigma2	WN variance
#
#	$aicc	Akaike informtion criterion corrected
#
#	$se.phi	0
#
#	$se.theta Standard errors of MA coefficients

ia = function(x,q,m=17) {

	if (q < 1)
		stop("q < 1")

	x = x - mean(x)

	# Innovations algorithm

	gamma = acvf(x,m+1)
	kappa = function(i,j) gamma[abs(i-j)+1]
	Theta = matrix(0,m,m)
	v = numeric(m+1)
	v[1] = kappa(1,1)
	for (n in 1:m) {
		for (k in 0:(n-1)) {
			if (k == 0)
				s = 0
			else
				s = sum(Theta[k,k:1] * Theta[n,n:(n-k+1)] * v[1:k])
			Theta[n,n-k] = (kappa(n+1,k+1) - s) / v[k+1]
		}
		s = sum(Theta[n,n:1]^2 * v[1:n])
		v[n+1] = kappa(n+1,n+1) - s
	}

	# MA coefficients

	theta = Theta[m,1:q]

	# Standard errors of MA coefficients

	n = length(x)
	se.theta = numeric(q)
	se.theta[1] = sqrt(1/n)
	if (q > 1)
		for (j in (2:q))
			se.theta[j] = sqrt((1+sum(Theta[m,1:(j-1)]^2))/n)

	# Estimated WN variance is not used

	sigma2 = v[m+1]

	a = list(
		phi=0,
		theta=theta,
		sigma2=NA,
		aicc=NA,
		se.phi=0,
		se.theta=se.theta)

	a = .innovation.update(x,a)

	return(a)
}

# Test the innovations algorithm (see p. 154 of ITSF)

.test.ia = function() {
	cat("ia ")
	book = c(
		1.9114, 1.1133, 0.4727, 0.6314, 0.5331, 0.6127, 0.4969, -0.0231,
		0.0568, -0.0064, 0.7594, -0.1757, 0.7666, 0.4801, -0.0792, -0.9563,
		0.2760)
	y = diff(dowj)
	a = ia(y,17)
	b = a$theta/a$se/1.96
	if (any(abs(b-book) > 1e-4))
		stop("fail")
	cat("ok\n")
}

# Hannan-Rissanen estimation algorithm for q > 0
#
# Arguments
#
#	x	Data
#
#	p	AR order
#
#	q	MA order
#
# Returns an ARMA model

hannan = function(x,p,q) {

	if (q < 1)
		stop("q < 1")

	n = length(x)
	x = x - mean(x)

	m = 20 + p + q #Reverse engineered from B&D ITSM
	k = max(p,q)   #Reverse engineered from B&D ITSM

	# Fit an AR(m)

	phi = ar.yw(x,aic=FALSE,order.max=m,demean=FALSE)$ar

	# Residuals of AR(m) process

	F1 = function(t) x[t] - sum(phi * x[(t-1):(t-m)])
	z = sapply((m+1):n,F1)
	z = c(numeric(m),z)

	# Design matrix (see p. 157 of B&D ITSF)

	F2 = function(i) z[(m+k+i-1):(m+k+i-q)]
	Z = sapply(1:(n-m-k),F2)
	Z = matrix(Z,nrow=n-m-k,ncol=q,byrow=TRUE)
	if (p > 0) {
		F3 = function(i) x[(m+k+i-1):(m+k+i-p)]
		X = sapply(1:(n-m-k),F3)
		X = matrix(X,nrow=n-m-k,ncol=p,byrow=TRUE)
		Z = cbind(X,Z)
	}

	# Regression

	G = qr.solve(qr(t(Z) %*% Z))
	b = G %*% t(Z) %*% x[(m+1+k):n]

	# Residuals

	xhat = Z %*% b
	e = x[(m+1+k):n] - xhat
	sigma2 = sum(e^2) / (n-m-k)
	se = sqrt(sigma2*diag(G))

	# AR coefficients

	if (p == 0) {
		phi = 0
		se.phi = 0
	} else {
		phi = b[1:p]
		se.phi = se[1:p]
	}

	# MA coefficients

	theta = b[p+(1:q)]
	se.theta = se[p+(1:q)]

	a = list(
		phi=phi,
		theta=theta,
		sigma2=NA,
		aicc=NA,
		se.phi=se.phi,
		se.theta=se.theta)

	a = .innovation.update(x,a)

	return(a)
}

# Estimate ARMA model parameters using maximum likelihood
#
# Arguments
#
#	x	Data
#
#	p	AR order
#
#	q	MA order
#
# Return value:
#
#	$phi	AR coefficients
#
#	$theta	MA coefficients
#
#	$sigma2	White noise variance
#
#	$aicc	Akaike information criterion corrected
#
#	$se.phi	Standard errors for phi
#
#	$se.theta Standard errors for theta

arma = function(x,p=0,q=0) {
	x = x - mean(x)
	w = getOption("warn")
	options(warn=-1) #Disable warnings
	a = arima(x,c(p,0,q))
	options(warn=w)
	c = as.vector(coef(a))
	v = as.vector(diag(vcov(a)))
	if (p == 0) {
		phi = 0
		se.phi = 0
	} else {
		phi = c[1:p]
		se.phi = v[1:p]
	}
	if (q == 0) {
		theta = 0
		se.theta = 0
	} else {
		theta = c[(p+1):(p+q)]
		se.theta = v[(p+1):(p+q)]
	}
	a = list(
		phi=phi,
		theta=theta,
		sigma2=NA,
		aicc=NA,
		se.phi=se.phi,
		se.theta=se.theta)
	a = .innovation.update(x,a)
	return(a)
}

# Returns the best ARMA model
#
# Arguments
#
#	x	Data
#
#	p	AR order sequence
#
#	q	MA order sequence

autofit = function(x,p=0:5,q=0:5) {
	a = list(aicc=Inf)
	for (j in p)
		for (k in q) {
			b = arma(x,j,k)
			if (b$aicc < a$aicc)
				a = b
		}
	return(a)
}

# Check an ARMA model for causality and invertibility
#
# Arguments
#
#	a	ARMA model (list of $phi, $theta)

check = function(a) {
	phi = c(1,-a$phi)
	theta = c(1,a$theta)
	if (all(abs(polyroot(phi)) > 1))
		cat("Causal\n")
	else
		cat("Non-Causal\n")
	if (all(abs(polyroot(theta)) > 1))
		cat("Invertible\n")
	else
		cat("Non-Invertible\n")
}

# Forecast future values of a time series
#
# Arguments
#
#	x	Observed data
#
#	xv	Transform vector (NULL for none)
#
#	a	ARMA model (list of $phi, $theta, $sigma2)
#
#	h	Number of predicted values
#
#	opt	Display option (0 = none, 1 = tabulate, 2 = plot and tabulate)
#
# Return value
#
#	$pred	Predicted values
#
#	$se	Standard errors (not included if log transform in xv)
#
#	$l	Lower 95% prediction bound
#
#	$u	Upper 95% prediction bound
#
# The transform vector xv is a vector of strings that specify a sequence of
# data transform functions.
#
# The transforms are applied left to right. Following this, the ARMA forecast
# is computed, then the inverse transforms are applied to the result.
#
# For example, the following transform takes the log of the data, then
# removes seasonality of period 12, then removes linear trend:
#
#	xv = c("log","season",12,"trend",1)

forecast = function(x,xv,a,h=10,opt=2) {

	f = .forecast.transform(x,xv,a,h,1)

	# Compute standard errors per (3.3.19) of ITSF

	if (!is.null(f$phi)) {
		psi = ma.inf(list(phi=f$phi,theta=a$theta),h)
		g = function(j) sum(psi[1:j]^2)
		se = sqrt(a$sigma2 * sapply(1:h,g))
		l = f$pred - 1.96*se
		u = f$pred + 1.96*se
		f = list(pred=f$pred,se=se,l=l,u=u)
	}

	if (opt > 0) {
		if (is.null(f$se))
			cat(" Step     Prediction    Lower Bound    Upper Bound\n")
		else
			cat(" Step     Prediction      sqrt(MSE)    Lower Bound    Upper Bound\n")
		for (i in 1:h) {
			cat(format(i,width=5))
			cat(format(f$pred[i],width=15))
			if (!is.null(f$se))
				cat(format(f$se[i],width=15))
			cat(format(f$l[i],width=15))
			cat(format(f$u[i],width=15))
			cat("\n")
		}
	}

	if (opt > 1)
		.plot.forecast(x,f)

	return(invisible(f))
}

# Transform the data, forecast, then invert the transform

.forecast.transform = function(x,xv,a,h,k) {

	if (k > length(xv))
		return(.forecast.arma(x,a,h))

	if (xv[k] == "diff")
		return(.forecast.diff(x,xv,a,h,k))

	if (xv[k] == "hr")
		return(.forecast.hr(x,xv,a,h,k))

	if (xv[k] == "log")
		return(.forecast.log(x,xv,a,h,k))

	if (xv[k] == "season")
		return(.forecast.season(x,xv,a,h,k))

	if (xv[k] == "trend")
		return(.forecast.trend(x,xv,a,h,k))

	stop(xv[k])
}

.forecast.arma = function(x,a,h) {

	n = length(x)

	mu = mean(x)
	x = x - mu

	dx = .innovations(x,a)
	dx = c(dx,numeric(h))
	x = c(x,numeric(h))

	phi = a$phi
	theta = a$theta

	p = length(phi)
	q = length(theta)

	# Forecast h steps ahead

	for (t in n:(n+h-1)) {
		A = sum(phi * x[t:(t-p+1)])
		B = sum(theta * dx[t:(t-q+1)])
		x[t+1] = A + B
	}

	pred = x[(n+1):(n+h)] + mu

	return(list(pred=pred,phi=phi))
}

# Difference the data, forecast, diffinv

.forecast.diff = function(x,xv,a,h,k) {

	n = length(x)

	lag = as.numeric(xv[k+1])

	y = diff(x,lag,1)
	f = .forecast.transform(y,xv,a,h,k+2)
	pred = diffinv(f$pred,lag,1,x[(n-lag+1):n])
	pred = pred[(lag+1):(lag+h)]

	# Update phi(z)

	# phi(z) = (1 - z^lag) * phi(z)

	if (is.null(f$phi))
		stop("diff before log")

	phi = f$phi
	phi = c(1,-phi)
	phi = c(phi,rep(0,lag)) - c(rep(0,lag),phi)
	phi = -phi[2:length(phi)]

	return(list(pred=pred,phi=phi))
}

# Subtract harmonic components, forecast, add them back

.forecast.hr = function(x,xv,a,h,k) {

	n = length(x)

	# Build the design matrix X

	n = length(x)
	X = rep(1,n+h)
	w = 0:(n-1+h)
	k = k+1
	while (k <= length(xv)) {
		if (regexpr("[^0-9]",xv[k]) > -1)
			break
		d = as.numeric(xv[k])
		X = cbind(X,cos(2*pi*w/d),sin(2*pi*w/d))
		k = k+1
	}

	# b is the regression coefficient vector

	b = solve.qr(qr(X[(1:n),]),x)

	# xhat is the harmonic component

	xhat = as.vector(X %*% b)

	# Subtract the harmonic component

	y = x - xhat[1:n]

	# Forecast

	f = .forecast.transform(y,xv,a,h,k)

	# Restore the harmonic component

	f$pred = f$pred + xhat[(n+1):(n+h)]

	if (is.null(f$phi)) {
		f$l = f$l + xhat[(n+1):(n+h)]
		f$u = f$u + xhat[(n+1):(n+h)]
	}

	return(f)
}

# Log, forecast, exponentiate

.forecast.log = function(x,xv,a,h,k) {

	f = .forecast.transform(log(x),xv,a,h,k+1)

	# Prediction bounds

	if (is.null(f$phi)) {
		l = f$l
		u = f$u
	} else {
		psi = ma.inf(list(phi=f$phi,theta=a$theta),h)
		g = function(j) sum(psi[1:j]^2)
		se = sqrt(a$sigma2 * sapply(1:h,g))
		l = f$pred - 1.96*se
		u = f$pred + 1.96*se
	}

	pred = exp(f$pred)
	l = exp(l)
	u = exp(u)

	return(list(pred=pred,l=l,u=u))
}

# Subtract seasonal component, forecast, add it back

.forecast.season = function(x,xv,a,h,k) {

	n = length(x)

	# d is the number of observations per season

	d = as.numeric(xv[k+1])

	# m is the estimated trend (ITSF p. 31)

	# Special case when d is even, see equation (1.5.12)

	q = floor(d/2)
	if (d == 2*q) {
		f = function(t) x[t-q]/2 + sum(x[(t-q+1):(t+q-1)]) + x[t+q]/2
		m = sapply((q+1):(n-q),f) / d
		m = c(rep(0,q),m,rep(0,q))
	} else
		m = smooth.ma(x,q)

	# w is the average deviation (ITSF p. 31)

	dx = x - m
	F = function(k) mean(dx[seq(k+q,n-q,d)])
	w = sapply(1:d,F)
	w = w - mean(w)

	# s is the seasonal component

	s = rep(w,len=n+q+h)[(q+1):(q+n+h)]

	# Subtract the seasonal component

	y = x - s[1:n]

	# Forecast

	f = .forecast.transform(y,xv,a,h,k+2)

	# Restore the seasonal component

	f$pred = f$pred + s[(n+1):(n+h)]

	if (is.null(f$phi)) {
		f$l = f$l + s[(n+1):(n+h)]
		f$u = f$u + s[(n+1):(n+h)]
	}

	return(f)
}

# Subtract trend component, forecast, add it back

.forecast.trend = function(x,xv,a,h,k) {

	n = length(x)

	# p is the order of the trend (1 linear, 2 quadratic, etc.)

	p = as.numeric(xv[k+1])

	# Build the design matrix X

	X = NULL
	for (j in (0:p))
		X = cbind(X,(1:(n+h))^j)

	# b is the regression coefficient vector

	b = solve.qr(qr(X[1:n,]),x)

	# xhat is the trend component

	xhat = as.vector(X %*% b)

	# Subtract the trend component

	y = x - xhat[1:n]

	# Forecast

	f = .forecast.transform(y,xv,a,h,k+2)

	# Restore the trend component

	f$pred = f$pred + xhat[(n+1):(n+h)]

	if (is.null(f$phi)) {
		f$l = f$l + xhat[(n+1):(n+h)]
		f$u = f$u + xhat[(n+1):(n+h)]
	}

	return(f)
}

# Estimates AR model parameters using the Yule-Walker method
#
# This is a wrapper for ar.yw()
#
# Arguments
#
#	x	Data
#
#	p	AR order
#
# Return value
#
#	$phi	AR coefficients
#
#	$theta	0
#
#	$sigma2	White noise variance
#
#	$aicc	Akaike information criterion corrected

yw = function(x,p) {
	x = x - mean(x)
	a = ar.yw(x,aic=FALSE,order.max=p)
	a = list(
		phi=a$ar,
		theta=0,
		sigma2=NA,
		aicc=NA,
		se.phi=diag(a$asy.var.coef),
		se.theta=0)
	a = .innovation.update(x,a)
	return(a)
}

# Estimates AR model parameters using the Burg method
#
# Arguments
#
#	x	Data
#
#	p	AR order
#
# Return value
#
#	$phi	AR coefficients
#
#	$theta	0
#
#	$sigma2	White noise variance
#
#	$aicc	Akaike information criterion corrected
#
#	$se.phi	Standard errors of AR coefficients
#
#	$se.theta 0

burg = function(x,p) {
	x = x - mean(x)
	a = ar.burg(x,aic=FALSE,order.max=p)
	a = list(
		phi=a$ar,
		theta=0,
		sigma2=NA,
		aicc=NA,
		se.phi=diag(a$asy.var.coef),
		se.theta=0)
	a = .innovation.update(x,a)
	return(a)
}

# Returns model residuals
#
# Arguments
#
#	x	Data
#
#	xv	Transform vector (NULL for none)
#
#	a	ARMA model (list of $phi, $theta or NULL for none)

Resid = function(x,xv=NULL,a=NULL) {
	y = x
	k = 1
	while (k < length(xv)) {
		if (xv[k] == "diff") {
			lag = as.numeric(xv[k+1])
			y = diff(y,lag,1)
			k = k+2
		} else if (xv[k] == "hr") {
			d = NULL
			k = k+1
			while (k <= length(xv)) {
				if (regexpr("[^0-9]",xv[k]) > -1)
					break
				d = c(d,as.numeric(xv[k]))
				k = k+1
			}
			y = y - hr(x,d)
		} else if (xv[k] == "log") {
			y = log(y)
			k = k+1
		} else if (xv[k] == "season") {
			d = as.numeric(xv[k+1])
			y = y - season(y,d)
			k = k+2
		} else if (xv[k] == "trend") {
			p = as.numeric(xv[k+1])
			y = y - trend(y,p)
			k = k+2
		} else
			stop(xv[k])
	}
	y = y - mean(y)
	if (!is.null(a)) {
		I = .innovation.kernel(y,a)
		y = (y - I$xhat) / sqrt(I$v)
	}
	return(y)
}

# Test for randomness of residuals
#
# Arguments
#
#	e	Residuals

test = function(e) {

	# Randomness tests

	cat("Null hypothesis: Residuals are iid noise.\n")
	cat("Test                        Distribution Statistic   p-value\n")
	.ljung_box_test(e)
	.mcleod_li_test(e)
	.turning_point_test(e)
	.difference_sign_test(e)
	.rank_test(e)

	# Plots

	n = length(e)
	h = min(40,n-1)
	op = par(mfrow=c(2,2))
	y = acf(e,lag.max=h,type="correlation",plot=FALSE)$acf
	plot(0:h,y,ylim=c(-1,1),type="h",xlab="Lag",ylab="",main="ACF",col="blue")
	n = length(e)
	b = 1.96/sqrt(n)
	abline(h=b,lty=2)
	abline(h=-b,lty=2)
	abline(h=0)
	y = acf(e,lag.max=h,type="partial",plot=FALSE)$acf
	y = c(1,y)
	plot(0:h,y,ylim=c(-1,1),type="h",xlab="Lag",ylab="",main="PACF",col="blue")
	abline(h=b,lty=2)
	abline(h=-b,lty=2)
	abline(h=0)
	plot(e,ylab="",main="Residuals")
	abline(h=0)
	qqnorm(e)
	par(op)
}

.ljung_box_test = function(y,h=20) {
	n = length(y)
	ybar = mean(y)
	acvf = function(h) sum((y[1:(n-h)]-ybar)*(y[(h+1):n]-ybar))/n
	gamma = sapply(0:h,acvf)
	rho = gamma/gamma[1]
	Q = n*(n+2)*sum(rho[2:(h+1)]^2/(n-(1:h)))
	pval = 1-pchisq(Q,h)
	.test.out(
		"Ljung-Box Q",
		Q,
		paste("Q ~ chisq(",h,")",sep=""),
		pval
	)
}

.mcleod_li_test = function(y,h=20) {
	y = y^2
	n = length(y)
	ybar = mean(y)
	acvf = function(h) sum((y[1:(n-h)]-ybar)*(y[(h+1):n]-ybar))/n
	gamma = sapply(0:h,acvf)
	rho = gamma/gamma[1]
	Q = n*(n+2)*sum(rho[2:(h+1)]^2/(n-(1:h)))
	pval = 1-pchisq(Q,h)
	.test.out(
		"McLeod-Li Q",
		Q,
		paste("Q ~ chisq(",h,")",sep=""),
		pval
	)
}

.turning_point_test = function(y) {
	n = length(y)
	a = y[1:(n-2)] < y[2:(n-1)] & y[2:(n-1)] > y[3:n]
	b = y[1:(n-2)] > y[2:(n-1)] & y[2:(n-1)] < y[3:n]
	T = sum(a)+sum(b)
	mu = 2*(n-2)/3
	sd = sqrt((16*n-29)/90)
	Z = (T-mu)/sd
	pval = 2*(1-pnorm(abs(Z),0,1))
	.test.out(
		"Turning points T",
		T,
		paste("(T-",round(mu,1),")/",round(sd,1)," ~ N(0,1)",sep=""),
		pval
	)
}

.difference_sign_test = function(y) {
	n = length(y)
	S = sum(y[2:n] > y[1:(n-1)])
	mu = (n-1)/2
	sd = sqrt((n+1)/12)
	Z = (S-mu)/sd
	pval = 2*(1-pnorm(abs(Z),0,1))
	.test.out(
		"Diff signs S",
		S,
		paste("(S-",round(mu,1),")/",round(sd,1)," ~ N(0,1)",sep=""),
		pval
	)
}

.rank_test = function(y) {
	n = length(y)
	f = function(j) sum(y[j:n] > y[j-1])
	P = sum(sapply(2:n,f))
	mu = n*(n-1)/4
	sd = sqrt(n*(n-1)*(2*n+5)/72)
	Z = (P-mu)/sd
	pval = 2*(1-pnorm(abs(Z),0,1))
	.test.out(
		"Rank P",
		P,
		paste("(P-",round(mu,1),")/",round(sd,1)," ~ N(0,1)",sep=""),
		pval
	)
}

.test.out = function(name,stat,dist,pval) {
	cat(name)
	cat(format(dist,width=40-nchar(name),justify="right"))
	cat(format(round(stat,2),width=10))
	cat(format(round(pval,4),width=10))
	if (pval < 0.05)
		cat(" *")
	cat("\n")
}

# Returns the innovations for an ARMA model
#
# Arguments
#
#	x	Data with mean subtracted off
#
#	a	ARMA model (list of $phi, $theta, $sigma2)

.innovations = function(x,a) {
	xhat = .innovation.kernel(x,a)$xhat
	return(x-xhat)
}

# Calculate xhat and v per ITSF
#
# Arguments
#
#	x	Data with mean subtracted off
#
#	a	ARMA model (list of $phi, $theta, $sigma2)
#
# Return value
#
#	$xhat	Estimate of data
#
#	$v	Mean square error

.innovation.kernel = function(x,a) {

	# Compute autocovariance kappa(i,j) per equation (3.3.3)

	# Optimized for i >= j and j > 0

	kappa = function(i,j) {
		if (j > m)
			return(sum(theta_r[1:(q+1)] * theta_r[(i-j+1):(i-j+q+1)]))
		else if (i > 2*m)
			return(0)
		else if (i > m)
			return((gamma[i-j+1] - sum(phi * gamma[abs(seq(1-i+j,p-i+j))+1]))/sigma2)
		else
			return(gamma[i-j+1]/sigma2)
	}

	phi = a$phi
	theta = a$theta
	sigma2 = a$sigma2

	N = length(x)

	theta_r = c(1,theta,numeric(N))

	# Autocovariance of the model

	gamma = aacvf(a,N-1)

	# Innovations algorithm

	p = length(phi)
	q = length(theta)
	m = max(p,q)

	Theta = matrix(0,N-1,N-1)
	v = numeric(N)
	v[1] = kappa(1,1)
	for (n in 1:(N-1)) {
		for (k in 0:(n-1)) {
			u = kappa(n+1,k+1)
			s = 0
			if (k > 0) s = sum(Theta[k,k:1] * Theta[n,n:(n-k+1)] * v[1:k])
			Theta[n,n-k] = (u - s) / v[k+1]
		}
		s = sum(Theta[n,n:1]^2 * v[1:n])
		v[n+1] = kappa(n+1,n+1) - s
	}

	# Compute xhat per equation (5.2.7)

	xhat = numeric(N)
	if (m > 1)
		for (n in 1:(m-1))
			xhat[n+1] = sum(Theta[n,1:n]*(x[n:1]-xhat[n:1]))
	for (n in m:(N-1)) {
		A = sum(phi*x[n:(n-p+1)])
		B = sum(Theta[n,1:q]*(x[n:(n-q+1)]-xhat[n:(n-q+1)]))
		xhat[n+1] = A + B
	}

	return(list(xhat=xhat,v=v))
}

# Update the AICC and WN variance of an ARMA model
#
# Arguments
#
#	x	Data with mean subtracted off
#
#	a	ARMA model (list of $phi, $theta)
#
# Return value
#
#	ARMA model with $aicc and $sigma2 added

.innovation.update = function(x,a) {

	a$sigma2 = 1
	I = .innovation.kernel(x,a)
	xhat = I$xhat
	v = I$v

	n = length(x)
	sigma2 = sum((x-xhat)^2/v)/n

	if (all(a$phi==0))
		p = 0
	else
		p = length(a$phi)

	if (all(a$theta==0))
		q = 0
	else
		q = length(a$theta)

	loglike = -(n/2)*log(2*pi*sigma2) - sum(log(v))/2 - n/2
	aicc = -2*loglike + 2*(p+q+1)*n/(n-p-q-2)

	a$sigma2 = sigma2
	a$aicc = aicc

	return(a)
}

# Forecast using ARAR algorithm
#
# Arguments
#
#	y	Data
#
#	h	Number of steps ahead
#
#	opt	Option (0 silent, 1 tabulate, 2 plot and tabulate)
#
# Return value
#
#	$pred	Predicted values
#
#	$se	Standard errors
#
#	$l	Lower bounds
#
#	$u	Upper bounds

arar = function(y,h=10,opt=2) {

	# Save y

	Y = y

	# Memory shortening step

	psi = 1

	f = function(tau) sum(y[(tau+1):n]*y[1:(n-tau)])/sum(y[1:(n-tau)]^2)
	g = function(tau) sum((y[(tau+1):n]-phi[tau]*y[1:(n-tau)])^2)/sum(y[(tau+1):n]^2)

	for (k in 1:3) {
		n = length(y)
		phi = sapply(1:15,f)
		err = sapply(1:15,g)
		tau = which.min(err)
		if (err[tau] <= 8/n | (phi[tau] >= 0.93 & tau > 2)) {
			y = y[(tau+1):n] - phi[tau]*y[1:(n-tau)]
			psi = c(psi,numeric(tau)) - phi[tau]*c(numeric(tau),psi)
		} else if (phi[tau] >= 0.93) {
			A = matrix(0,2,2)
			A[1,1] = sum(y[2:(n-1)]^2)
			A[1,2] = sum(y[1:(n-2)]*y[2:(n-1)])
			A[2,1] = sum(y[1:(n-2)]*y[2:(n-1)])
			A[2,2] = sum(y[1:(n-2)]^2)
			b = c(sum(y[3:n]*y[2:(n-1)]),sum(y[3:n]*y[1:(n-2)]))
			phi = qr.solve(qr(A),b)
			y = y[3:n] - phi[1]*y[2:(n-1)] - phi[2]*y[1:(n-2)]
			psi = c(psi,0,0) - phi[1]*c(0,psi,0) - phi[2]*c(0,0,psi)
		} else
			break()
	}

	S = y
	Sbar = mean(S)
	X = S - Sbar
	gamma = acvf(X)

	# Restore y

	y = Y

	# Find best lags

	A = matrix(gamma[1],4,4)
	b = numeric(4)
	best.sigma2 = Inf
	m = 26
	for (i in 2:(m-2)) for (j in (i+1):(m-1)) for (k in (j+1):m) {
		A[1,2] = A[2,1] = gamma[i]
		A[1,3] = A[3,1] = gamma[j]
		A[2,3] = A[3,2] = gamma[j-i+1]
		A[1,4] = A[4,1] = gamma[k]
		A[2,4] = A[4,2] = gamma[k-i+1]
		A[3,4] = A[4,3] = gamma[k-j+1]
		b[1] = gamma[2]
		b[2] = gamma[i+1]
		b[3] = gamma[j+1]
		b[4] = gamma[k+1]
		phi = qr.solve(qr(A),b)
		sigma2 = gamma[1]-phi[1]*gamma[2]-phi[2]*gamma[i+1]-phi[3]*gamma[j+1]-phi[4]*gamma[k+1]
		if (sigma2 < best.sigma2) {
			best.sigma2 = sigma2
			best.phi = phi
			best.lag = c(1,i,j,k)
		}
	}

	i = best.lag[2]
	j = best.lag[3]
	k = best.lag[4]

	phi = best.phi
	sigma2 = best.sigma2

	xi = c(psi,numeric(k))-
		phi[1]*c(0,psi,numeric(k-1))-
		phi[2]*c(numeric(i),psi,numeric(k-i))-
		phi[3]*c(numeric(j),psi,numeric(k-j))-
		phi[4]*c(numeric(k),psi)

	# Forecast

	n = length(y)
	k = length(xi)
	y = c(y,numeric(h))
	c = (1-sum(phi))*Sbar
	for (i in 1:h)
		y[n+i] = -sum(xi[2:k]*y[n+i+1-(2:k)]) + c
	pred = y[n+(1:h)]
	y = y[1:n]

	# Extend filter coefficients if necessary

	if (h > k)
		xi = c(xi,numeric(h-k))

	# See equation (9.1.9)

	tau = numeric(h)
	tau[1] = 1
	if (h > 1)
		for (j in 1:(h-1))
			tau[j+1] = -sum(tau[1:j]*xi[(j:1)+1])

	# Standard errors from (9.1.8), except correct typo

	F = function(j) sum(tau[1:j]^2)
	se = sqrt(sigma2*sapply(1:h,F))

	f = list(pred=pred,se=se,l=pred-1.96*se,u=pred+1.96*se)

	if (opt > 0) {
		cat("Optimal lags",best.lag,"\n")
		cat("Optimal coeffs",best.phi,"\n")
		cat("WN Variance",best.sigma2,"\n")
		cat("Filter", xi,"\n\n")
		cat(" Step     Prediction      sqrt(MSE)    Lower Bound    Upper Bound\n")
		for (i in 1:h) {
			cat(format(i,width=5))
			cat(format(f$pred[i],width=15))
			cat(format(f$se[i],width=15))
			cat(format(f$l[i],width=15))
			cat(format(f$u[i],width=15))
			cat("\n")
		}
	}

	if (opt > 1)
		.plot.forecast(y,f)

	return(invisible(f))
}

.test.arar.deaths = function() {
	cat("test arar (deaths) ")
	f = arar(deaths,opt=0)
	pred = c(8167.8,7195.8,7982.0,8283.5,9144.1,9464.9,10541.0,9640.8,8902.7,9096.7)
	se = c(323.35,375.68,392.34,414.78,431.69,441.90,449.82,455.98,460.45,463.77)
	l = c(7534.1,6459.4,7213.1,7470.6,8298.0,8598.8,9659.7,8747.1,8000.3,8187.7)
	u = c(8801.6,7932.1,8751.0,9096.5,9990.2,10331.0,11423.0,10535.0,9805.2,10006.0)
	.test.arar.kernel(f,pred,se,l,u,0.0001)
	cat("ok\n")
}

# trings data has medium memory

.test.arar.trings = function() {
	cat("test arar (trings) ")
	f = arar(.trings,opt=0)
	pred = c(0.97008,0.99522,1.04256,0.99059,1.02546,1.02129,1.01924,1.00636,1.01629,1.01310)
	se = c(0.18340,0.20841,0.22528,0.23397,0.24923,0.25677,0.26742,0.27785,0.28886,0.29755)
	l = c(0.61061,0.58675,0.60102,0.53203,0.53697,0.51802,0.49511,0.46178,0.45014,0.42992)
	u = c(1.32954,1.40369,1.48410,1.44916,1.51394,1.52455,1.54337,1.55095,1.58244,1.59628)
	.test.arar.kernel(f,pred,se,l,u,0.0001)
	cat("ok\n")
}

# For arar, ITSM-R sigma2 is different from B&D ITSM so divide by sigma2 to compare.

.test.arar.kernel = function(f,pred,se,l,u,eps) {
	e = (f$pred - pred) / pred
	if (any(abs(e) > eps)) {
		cat("\n")
		print(f$pred)
		print(pred)
		print(e)
		stop("fail pred\n")
	}
	e = (f$l - (f$pred - 1.96*f$se)) / f$l
	if (any(abs(e) > eps)) {
		cat("\n")
		print(f$l)
		print(e)
		stop("fail lower bound\n")
	}
	e = (f$u - (f$pred + 1.96*f$se)) / f$u
	if (any(abs(e) > eps)) {
		cat("\n")
		print(f$u)
		print(e)
		stop("fail upper bound\n")
	}
	f$se = f$se/f$se[1]
	se = se/se[1]
	e = (f$se - se) / se
	if (any(abs(e) > eps)) {
		cat("\n")
		print(f$se)
		print(se)
		print(e)
		stop("fail se\n")
	}
}

.test.forecast.lake = function() {
	cat("test forecast (lake) ")
	xv = c("diff",1)
	e = Resid(lake,xv)
	a = arma(e,0,0)
	f = forecast(lake,xv,a,opt=0)
	pred = c(9.95567,9.95134,9.94701,9.94268,9.93835,9.93402,9.92969,9.92536,9.92103,9.9167)
	se = c(0.74518,1.05384,1.29069,1.49036,1.66627,1.82531,1.97156,2.10768,2.23554,2.35646)
	l = c(8.49515,7.88585,7.41731,7.02163,6.67252,6.35648,6.06551,5.79438,5.53946,5.29812)
	u = c(11.41619,12.01683,12.47671,12.86373,13.20418,13.51156,13.79387,14.05634,14.3026,14.53528)
	.test.forecast.kernel(f,pred,se,l,u,0.0001)
	cat("ok\n")
}

.test.forecast.dowj = function() {
	cat("test forecast (dowj) ")
	xv = c("diff",1)
	e = Resid(dowj,xv)
	a = arma(e,1,0)
	f = forecast(dowj,xv,a,opt=0)
	pred = c(120.96,120.91,120.97,121.06,121.18,121.31,121.44,121.57,121.7,121.84)
	se = c(0.38146,0.67099,0.91921,1.13299,1.32016,1.48703,1.63825,1.77718,1.90622,2.02716)
	l = c(120.21,119.6,119.16,118.84,118.59,118.39,118.23,118.09,117.97,117.86)
	u = c(121.71,122.23,122.77,123.28,123.77,124.22,124.65,125.05,125.44,125.81)
	.test.forecast.kernel(f,pred,se,l,u,0.005)
	cat("ok\n")
}

.test.forecast.wine = function() {
	cat("test forecast (wine) ")
	xv = c("log","diff",12)
	e = Resid(wine,xv)
	a = arma(e,1,0)
	f = forecast(wine,xv,a,opt=0)
	pred = c(2323.8,2456.5,1079.4,1783.2,1758,1632.6,1967.6,2025.4,3125.9,2753.4)
	l = c(1780.9,1853.9,813.19,1343.1,1324.1,1229.7,1482,1525.5,2354.4,2073.8)
	u = c(3032.3,3254.9,1432.8,2367.5,2334.1,2167.6,2612.4,2689.1,4150.2,3655.6)
	.test.forecast.kernel(f,pred,NULL,l,u,0.0001)
	cat("ok\n")
}

.test.forecast.deaths = function() {
	cat("test forecast (deaths) ")
	xv = c("diff",12)
	e = Resid(deaths,xv)
	a = arma(e,1,1)
	f = forecast(deaths,xv,a,opt=0)
	pred = c(8175.4,7193.9,8058.2,8364,9320.2,9611.6,10636,9955.2,9216.3,9155.9)
	se = c(351.83,394.64,428,454.69,476.41,494.3,509.16,521.57,532,540.8)
	l = c(7485.8,6420.4,7219.3,7472.8,8386.4,8642.7,9638,8933,8173.6,8096)
	u = c(8864.9,7967.4,8897,9255.2,10254,10580,11634,10978,10259,10216)
	.test.forecast.kernel(f,pred,NULL,l,u,0.001)
	cat("ok\n")
}

.test.forecast.sunspots = function() {
	cat("test forecast (sunspots) ")
	a = arma(Sunspots,2,1)
	f = forecast(Sunspots,NULL,a,opt=0)
	pred = c(88.3043,82.43004,67.2167,51.87652,41.61624,37.64887,38.54133,41.85841,45.42085,47.92461)
	se = c(14.62727,27.71399,34.55017,36.59889,36.74043,36.84864,37.22249,37.53885,37.65907,37.67121)
	l = c(59.63538,28.11162,-0.50039,-19.85599,-30.39368,-34.57313,-34.4134,-31.71638,-28.38957,-25.9096)
	u = c(116.97,136.75,134.93,123.61,113.63,109.87,111.5,115.43,119.23,121.76)
	.test.forecast.kernel(f,pred,se,l,u,0.05) #To accomodate step 3 lower bound
	cat("ok\n")
}

.test.forecast.kernel = function(f,pred,se,l,u,eps) {
	e = (f$pred - pred) / pred
	if (any(abs(e) > eps)) {
		cat("\n")
		print(f$pred)
		print(pred)
		print(e)
		stop("fail pred\n")
	}
	if (!is.null(f$se)) {
		e = (f$se - se) / se
		if (any(abs(e) > eps)) {
			cat("\n")
			print(f$se)
			print(se)
			print(e)
			stop("fail se\n")
		}
	}
	e = (f$l - l) / l
	if (any(abs(e) > eps)) {
		cat("\n")
		print(f$l)
		print(l)
		print(e)
		stop("fail lower bound\n")
	}
	e = (f$u - u) / u
	if (any(abs(e) > eps)) {
		cat("\n")
		print(f$u)
		print(u)
		print(e)
		stop("fail upper bound\n")
	}
}

# Compute the periodogram
#
# Argument
#
#	x	Data
#
#	q	Moving average filter order (can be a vector)
#
#	opt	Option
#			= 0 Silent
#			= 1 Plot periodogram only
#			= 2 Plot both periodogram and filter
#
# Returns a vector of estimated frequency amplitudes
#
#	f(lambda)
#
# where
#
#	lambda = pi * (1:m)/m
#
#	m = floor(length(x)/2)

periodogram = function(x,q=0,opt=2) {
	n = length(x)
	t = 1:n
	m = floor(n/2)
	k = pi*(1:m)/m
	I = function(lambda) (sum(x*cos(t*lambda))^2+sum(x*sin(t*lambda))^2)/n
	y = sapply(k,I) / (2*pi)
	w = 1
	if (any(q > 0)) {
		for (j in q)
			if (j > 0)
				w = .fcomp(w,j)
		j = floor(length(w)/2)
		y = c(rep(y[1],j),y,rep(y[m],j))
		j = 0:(2*j)
		f = function(t) sum(y[t+j] * w)
		y = sapply(1:m,f)
	}
	if (opt > 0)
		if (all(q == 0)) {
			main = "Periodogram/2pi"
			plot(k,y,type="o",col="blue",xlab="",ylab="",main=main)
		} else {
			cat("Filter ",w,"\n")
			main = "(Smoothed Periodogram)/2pi"
			if (opt == 1)
				plot(k,y,type="o",col="blue",xlab="",ylab="",main=main)
			else {
				op = par(mfrow=c(2,1))
				plot(k,y,type="o",col="blue",xlab="",ylab="",main=main)
				plot(w,type="h",col="blue",xlab="",ylab="",main="Smoothing Filter")
				par(op)
			}
		}
	return(invisible(y))
}

# Moving average filter composition
#
#	w	Weights of current filter
#
#	q	Order of next MA filter
#
# Returns a vector of the new filter weights

.fcomp = function(w,q) {
	n = length(w)
	w = c(rep(0,2*q),w,rep(0,2*q))
	k = 0:(2*q)
	f = function(t) sum(w[t+k])
	y = sapply(1:(2*q+n),f) / (2*q+1)
	return(y)
}

.test.periodogram = function() {
	cat("test periodogram ")
	y = periodogram(Sunspots,opt=0)
	yy = c(
		1346.5,1382.8,40.708,174.97,46.414,56.393,959.75,
		710.34,1514.5,2174.6,189.83,1202.1,127.33,221.34,16.062,
		36.825,188.78,125.13,49.743,24.446,86.523,69.833,46.737,
		52.925,1.364,41.845,15.408,4.5947,8.0343,0.76497,10.794,
		6.2038,8.1066,21.795,0.16017,5.225,7.1887,0.12595,0.25181,
		1.2706,1.7537,0.87732,1.1386,0.65564,2.407,3.1375,5.7411,
		2.2302,1.0778,4.8144)
	e = (y - yy)/yy
	if (any(abs(e) > 1e-4)) {
		print(e)
		stop("fail 1st test")
	}
	y = periodogram(Sunspots,1,opt=0)
	yy = c(
		1358.571,923.3199,532.8227,87.36368,92.59193,354.1869,
		575.4954,1061.524,1466.488,1292.985,1188.845,506.4071,
		516.9095,121.5766,91.4083,80.55662,116.9131,121.2190,
		66.43994,53.5705,60.26729,67.69768,56.49815,33.67515,
		32.04469,19.53912,20.61605,9.34567,4.464658,6.5312,5.921034,
		8.368245,12.03525,10.02073,9.060195,4.191301,4.179895,
		2.522161,0.5494647,1.092052,1.300555,1.256559,0.8905358,
		1.400438,2.066728,3.76187,3.702941,3.016359,2.707482,
		3.568881)
	e = (y - yy)/yy
	if (any(abs(e) > 1e-4)) {
		print(e)
		stop("fail 2nd test")
	}
	y = periodogram(Sunspots,2,opt=0)
	yy = c(
		1092.576,858.2777,598.2685,340.255,255.6475,389.5738,
		657.4761,1083.122,1109.810,1158.271,1041.669,783.0407,
		351.3243,320.7231,118.0675,117.6279,83.30881,84.9855,
		94.9252,71.13521,55.45624,56.09268,51.47634,42.54076,
		31.65577,23.22737,14.24927,14.12948,7.919271,6.078421,
		6.780792,9.533022,9.412061,8.298189,8.495183,6.899051,
		2.590333,2.812426,2.118165,0.8558852,1.058425,1.139190,
		1.366469,1.64323,2.615979,2.834298,2.918724,3.400206,
		3.73559,3.550264)
	e = (y - yy)/yy
	if (any(abs(e) > 1e-4)){
		print(e)
		stop("fail 3rd test")
	}
	y = periodogram(Sunspots,c(1,2),opt=0)
	yy = c(
		1101.527,849.7074,598.9337,398.057,328.4921,434.2325,
		710.0573,950.136,1117.068,1103.25,994.327,725.3447,485.0294,
		263.3716,185.4728,106.3347,95.3074,87.73983,83.68197,
		73.83888,60.89471,54.34175,50.03659,41.89096,32.47463,
		23.04414,17.20204,12.09934,9.375723,6.926161,7.464078,
		8.575292,9.08109,8.735144,7.897474,5.994856,4.100603,
		2.506975,1.928826,1.344158,1.017833,1.188028,1.382963,
		1.875226,2.364502,2.789667,3.051076,3.351507,3.56202,
		3.784319)
	e = (y - yy)/yy
	if (any(abs(e) > 1e-4)) {
		print(e)
		stop("fail 4th test")
	}
	cat("ok\n")
}
