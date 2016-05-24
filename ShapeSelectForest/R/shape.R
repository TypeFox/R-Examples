getedf0 = function(x, flat = TRUE, dec = TRUE, jp = TRUE, invee = TRUE, vee = TRUE, inc = TRUE, db = TRUE, nsim = 1e+3, random = FALSE, msg = TRUE) {
	n = length(x)
	k0 = 4 + round(n^(1 / 7)) + 2
	k1 = k0
	k = NULL
	while (is.null(k)) { 
		pts = floor((n - 2) / (k1 - 1))
		rem_pts = (n - 2) %% (k1 - 1)
		if (pts > 2) {
			k = k1
		} else if (pts == 2) {
			if (rem_pts / (k1 - 1) >= 1) {
				k = k1
			} else {
				k1 = k1 - 1	
			}
		} else {
			k1 = k1 - 1
		}
	}
	if (nsim < 1e+3) {
		stop ("We need at least 1000 simulations to get edf0!")
	}
	nloop = nsim
	xs = sort(x)
	x = (xs - min(xs)) / (max(xs) - min(xs))
	kobs = 1:k

	ans = bqspl(x, k, knots = NULL, pic = FALSE)
	knots = ans$knots
	slopes = ans$slopes
	delta = ans$bmat

	#qv = t(delta) %*% delta
	qv = crossprod(delta)	
	m0 = length(delta) / n	
	#umat0 = chol(t(delta) %*% delta)
	umat0 = chol(qv)
	uinv0 = solve(umat0)
	pmult0 = t(uinv0) %*% t(delta)
## constraint matrix for decreasing
	amat1 = -slopes %*% uinv0
if (random) {
	sz = sample(1:6, size = 1)
	loc = sample(1:7, size = sz, replace = FALSE) 
	bool = rep(TRUE, 7)
	bool[loc] = FALSE
	flat = bool[1]; dec = bool[2]; jp = bool[3]; invee = bool[4]; vee = bool[5]; inc = bool[6]; db = bool[7]
} else {
	bool = c(flat, dec, jp, invee, vee, inc, db)
}
edf0 = NULL
if (flat) {
	edf0 = c(edf0, 1)
}
if (dec) {
## get decreasing expected DF0 
	sd = 0
	for (iloop in 1:nloop){
		ys = rnorm(n)
		z = pmult0 %*% ys
		ans = coneA(z, amat1, msg = msg)
		sd = sd + ans$df
	}
	edf2 = sd / nloop
	edf0 = c(edf0, edf2)
}
if (jp) {
## get jump expected DF0 
	mj = length(delta) / n + 2
	sd = 0
	for (iloop in 1:nloop) {
		ys = rnorm(n)
		minsse = sum((ys - mean(ys))^2)
		for (i in 1:(n - 1)) {
			smat = matrix(0, nrow = k + 3, ncol = mj)
			smat[1:k, 1:(mj - 2)] = -slopes
			smat[k + 1, mj - 1] = 1

			djump = matrix(0, nrow = n, ncol = mj)
			djump[ ,1:(mj - 2)] = delta

			djump[1:i, mj - 1] = (2 * i + 1) / 2 - n; djump[(i + 1):n, mj - 1] = (2 * i + 1) / 2
			djump[1:i, mj] = 0; djump[(i + 1):n, mj] = x[x > (x[i] + x[i+1]) / 2] - (x[i] + x[i+1]) / 2
			djump[, mj] = djump[, mj] - mean(djump[, mj])

			kn = 1:(k + 3) < 0
			kn[1:k] = knots > (x[i] + x[i+1])/2

			smat[, mj] = 0; smat[kn, mj] = -1; smat[k + 2, mj] = -1
			ansi = -sl((x[i] + x[i+1]) / 2, knots, slopes)	
			smat[k + 2, 1:m0] = ansi 
			smat[k + 3, 1:m0] = ansi  	

			use = 1:(k + 3) > 0
			kp = min(knots[knots > (x[i] + x[i+1]) / 2])
			if (sum(x > (x[i] + x[i+1]) / 2 & x <= kp) == 0) {use[k + 2] = FALSE}
			kp = max(knots[knots < (x[i] + x[i+1]) / 2])
			if (sum(x < (x[i] + x[i+1]) / 2 & x >= kp) == 0) {use[k + 3] = FALSE}
			smat = smat[use, ]

			if (i == 1 | i == (n - 1)) {
				smat = smat[, -mj, drop = FALSE]
				djump = djump[, -mj, drop = FALSE] 
			}

			#umat = chol(t(djump) %*% djump)
			umat = chol(crossprod(djump))			
			uinv = solve(umat)
			pmult = t(uinv) %*% t(djump)
			z = pmult %*% ys
			amat = smat %*% uinv
			fiti = coneA(z, amat, msg = msg)
			phat = fiti$thetahat
			theta = djump %*% uinv %*% phat
			ssei = sum((ys - theta)^2)
			if (ssei < minsse) {ijump = i; sse3 = ssei; thb = theta; df = fiti$df; minsse = ssei}
		}
		sd = sd + df
#		print(100*iloop+df)
	}
	edf3 = sd / nloop + 1  # add one for jump point parameter
	edf0 = c(edf0, edf3)
}
if (vee | invee) {
## get vee and inverted vee expected DF0 	
	sd = 0
	for (iloop in 1:nloop) {
		ys = rnorm(n)
		#cv = t(delta) %*% ys
		cv = crossprod(delta, ys)
		minsse = sum((ys - mean(ys))^2)
		av = slopes
		for (j in 2:k) {
			av1 = -av
			av1[j:k,] = -av1[j:k,]
			qans = qprog(qv, cv, av1, 1:k*0, msg = msg)
			theta = delta %*% qans$thetahat
			sse = sum((ys - theta)^2)
			if (sse < minsse) {cch = j; thb = theta; minsse = sse; df = qans$df}
		}
		sd = sd + df
	}
	edf4 = sd / nloop + 1
	nv = sum(c(vee, invee))
	edf0 = c(edf0, rep(edf4, nv))
}
if (inc) {
	edf0 = c(edf0, 1.5)
}
#new: db == TRUE or FALSE
if (db) {
## get two jump expected DF0
	m0 = length(delta) / n
	mj =  m0 + 4
	sd = 0
	for (iloop in 1:nloop) {
		ys = rnorm(n)
		minsse = sum((ys - mean(ys))^2)
		rm_vec = NULL
		for (i in 1:(n - 4)) {
			for (j in (i + 3):(n - 1)) {
				smat = matrix(0, nrow = k + 6, ncol = mj)
				smat[1:k, 1:(mj - 4)] = -slopes
				smat[k + 1, mj - 3] = 1
				#smat[k + 4, mj - 1] = 1

				djump = matrix(0, nrow = n, ncol = mj)
				djump[ ,1:(mj - 4)] = delta

				djump[1:i, mj - 3] = (2 * i + 1) / 2 - n; djump[(i + 1):n, mj - 3] = (2 * i + 1) / 2
				djump[1:i, mj - 2] = 0; djump[(i + 1):n, mj - 2] = x[x > (x[i] + x[i + 1]) / 2] - (x[i] + x[i + 1]) / 2
				djump[, mj - 2] = djump[, mj - 2] - mean(djump[, mj - 2])
		
				djump[1:j, mj - 1] = (2 * j + 1) / 2 - n; djump[(j + 1):n, mj - 1] = (2 * j + 1) / 2
				djump[1:j, mj] = 0; djump[(j + 1):n, mj] = x[x > (x[j] + x[j + 1]) / 2] - (x[j] + x[j + 1]) / 2
				djump[, mj] = djump[, mj] - mean(djump[, mj])
			
				kn = 1:(k + 6) < 0
				kn[1:k] = knots > (x[i] + x[i+1]) / 2
				smat[kn, mj - 2] = -1 
				smat[k + 2, mj - 2] = -1
				ansi = -sl((x[i] + x[i+1]) / 2, knots, slopes)
				smat[k + 2, 1:m0] = ansi 
				smat[k + 3, 1:m0] = ansi 
			
				kn = 1:(k + 6) < 0
				kn[1:k] = knots > (x[j] + x[j+1]) / 2
				smat[kn, mj] = -1
				smat[k + 5, mj] = -1
				smat[k + 5, mj - 2] = -1
				ansj = -sl((x[j] + x[j+1]) / 2, knots, slopes)
				smat[k + 5, 1:m0] = ansj
				smat[k + 6, 1:m0] = ansj 
				smat[k + 6, mj - 2] = -1

				use = 1:(k + 6) > 0
				kp = min(knots[knots > (x[i] + x[i+1]) / 2])
				if (sum(x > (x[i] + x[i+1]) / 2 & x <= kp) == 0){use[k + 2] = FALSE}
				kp = max(knots[knots < (x[i] + x[i+1]) / 2])
				if (sum(x < (x[i] + x[i+1]) / 2 & x >= kp) == 0){use[k + 3] = FALSE}
				kp = min(knots[knots > (x[j] + x[j+1]) / 2])
				if (sum(x > (x[j] + x[j+1]) / 2 & x <= kp) == 0){use[k + 5] = FALSE}
				kp = max(knots[knots < (x[j] + x[j+1]) / 2])	
				if (sum(x < (x[j] + x[j+1]) / 2 & x >= kp) == 0){use[k + 6] = FALSE}
				smat = smat[use, ]
	
				if (i == 1) {
					rm_vec = c(rm_vec, mj - 2)
				}

				if (j == (n - 1)) {
					rm_vec = c(rm_vec, mj)
				}

				if (!is.null(rm_vec)) {
					smat = smat[, -rm_vec, drop = FALSE]
					djump = djump[, -rm_vec, drop = FALSE]
					rm_vec = NULL
				}

				#umat = chol(t(djump) %*% djump)
				umat = chol(crossprod(djump))
				uinv = solve(umat)
				pmult = t(uinv) %*% t(djump)
				z = pmult %*% ys
				amat = smat %*% uinv			
				fiti = coneA(z, amat, msg = msg)
				phat = fiti$thetahat
				theta = djump %*% uinv %*% phat
				ssei = sum((ys - theta)^2)
			
			if (ssei < minsse) {ijump = i; jjump = j; sse7 = ssei; thb = theta; df = fiti$df; minsse = ssei}
			}
		}
		sd = sd + df
	}
	edf5 = sd / nloop + 2
	edf0 = c(edf0, edf5)
}
## define edf0 for flat, decr, jump, inv, vee, incr, and 2jump respectively
#new: db == TRUE or FALSE
	#if (db) {
	#	ne = length(edf0)
	#	edf0 = c(edf0[1:(ne - 1)], 1.5, edf5)
		#edf0 = c(1, edf2, edf3, edf4, edf4, 1.5, edf5)
	#} else {
		#edf0 = c(1, edf2, edf3, edf4, edf4, 1.5)
	#	edf0 = c(edf0, 1.5)
	#}
	edf0	
}


## Main code
shape = function(x, ymat, infocrit = "CIC", flat = TRUE, dec = TRUE, jp = TRUE, invee = TRUE, vee = TRUE, inc = TRUE, db = TRUE, nsim = 1e+3, edf0 = NULL, get.edf0 = FALSE, random = FALSE, msg = TRUE) {
	data("edf0s", package = "ShapeSelectForest", envir = environment())
	edf0s = get("edf0s")
	ymat = as.matrix(ymat)
	x0 = x	
	xs = sort(x)
	n = length(xs)
	x = (xs - min(xs)) / (max(xs) - min(xs))
#new: any shape could be excluded
if (random) {
	sz = sample(1:6, size = 1)
	loc = sample(1:7, size = sz, replace = FALSE) 
	bool = rep(TRUE, 7)
	bool[loc] = FALSE
	flat = bool[1]; dec = bool[2]; jp = bool[3]; invee = bool[4]; vee = bool[5]; inc = bool[6]; db = bool[7]
} else {
	bool = c(flat, dec, jp, invee, vee, inc, db)
}
	shape = (1:7)[bool]
	nsh = sum(bool)
#new: prelim check with edf0
	#if (!is.null(edf0)) {
	#	if (length(edf0) < nsh) {
	#		stop ("Edf0 vector is wrong! Check to make sure that each shape has an edf0!")
	#	}
	#}
#new: check if equally spaced
	dist = round(diff(x), 10)
	if (all(dist == dist[1]) & n <= 40 & n >= 20) {
		if (is.null(edf0)) {
			#if (db) {
			#	edf0 = edf0s[(n - 20 + 1), ]
			#} else {
			#	edf0 = edf0s[(n - 20 + 1), (1:6)]
			#}
			edf0 = edf0s[(n - 20 + 1), bool] 
		} else {
			#if (db) {
			#	if (length(edf0) != 7) {
			#		stop ("Double-jump is allowed! The edf0 vector should have 7 elements!")
			#	}
			#} else {
			#	if (length(edf0) >= 6) {
			#		edf0 = edf0[1:6]
			#	} else {
			#		stop ("There should be >= 6 elements in the edf0 vector!")
			#	}
			#}
			if (length(edf0) > nsh) {
				edf0 = edf0[1:nsh]
			}
			if (length(edf0) < nsh) {			
				stop ("Edf0 vector is wrong! Check to make sure that each shape has an edf0!")
			}
		}
	} else {
		if (get.edf0) {
			if (n > 40) {
				print ("The x vector has more than 40 elements! The edf0 vector is being calculated!")
			} else if (n < 20) {
				print ("The x vector has less than 20 elements! The edf0 vector is being calculated!")
			} else if (!all(dist == dist[1])) {
				print ("The x vector is not equally spaced! The edf0 vector is being calculated!") 
			}
			edf0 = getedf0(x, flat = flat, dec = dec, jp = jp, invee = invee, vee = vee, inc = inc, db = db, nsim = nsim) 
		} else {
			if (is.null(edf0)) {
				if (n > 40) {
					stop ("The x vector has more than 40 elements! A edf0 vector should be provided or get.edf0 should be TRUE!")
				} else if (n < 20) {
					stop ("The x vector has less than 20 elements! A edf0 vector should be provided or get.edf0 should be TRUE!")
				} else if (!all(dist == dist[1])) {
					stop ("The x vector is not equally spaced! A edf0 vector should be provided or get.edf0 should be TRUE!")
				}
			} else {
				#if (db) {
				#	if (length(edf0) != 7) {
				#		stop ("Double-jump is allowed! The edf0 vector should have 7 elements!")
				#	}
				#} else {
				#	if (length(edf0) >= 6) {
				#		edf0 = edf0[1:6]
				#	} else {
				#		stop ("There should be >= 6 elements in the edf0 vector!")
				#	}
				#}
				if (length(edf0) > nsh) {
					edf0 = edf0[1:nsh]
				}
				if (length(edf0) < nsh) {					
					stop ("Edf0 vector is wrong! Check to make sure that each shape has an edf0!")
				}
			}
		}
	}
	ny = length(ymat) / n
	#k = 4 + round(n^(1 / 5)) + 2
	k0 = 4 + round(n^(1 / 7)) + 2
	k1 = k0
	k = NULL
	while (is.null(k)) { 
		pts = floor((n - 2) / (k1 - 1))
		rem_pts = (n - 2) %% (k1 - 1)
		if (pts > 2) {
			k = k1
		} else if (pts == 2) {
			if (rem_pts / (k1 - 1) >= 1) {
				k = k1
			} else {
				k1 = k1 - 1	
			}
		} else {
			k1 = k1 - 1
		}
	}

	kobs = 1:k
	ans = bqspl(x, k, knots = NULL, pic = FALSE)
	knots = ans$knots
	slopes = ans$slopes
	delta = ans$bmat

	qv = crossprod(delta)
	m0 = length(delta) / n
	#umat0 = chol(t(delta) %*% delta)
	umat0 = chol(qv)	
	uinv0 = solve(umat0)
	pmult0 = t(uinv0) %*% t(delta)
## constraint matrix for decreasing
	amat1 = -slopes %*% uinv0
## loop through ymat to choose shapes
#new: include a double-jump or not
	#if (db) {
	#	nsh = 7
	#} else {
	#	nsh = 6
	#}
	ic = matrix(nrow = ny, ncol = nsh)
	thetab = matrix(nrow = ny, ncol = n)
	fit = matrix(nrow = nsh, ncol = n)
	shb = 1:ny
	#shape = 1:nsh; shb = 1:ny
#new:
	ijps0 = jjps0 = 1:nsh*0
	m_is0 = m_js0 = 1:nsh*0
	ijps = jjps = list()
	m_is = m_js = list()
	bs = bhat = list()
	for (s in 1:ny) {
ish = 0 
		y = ymat[, s]
sse1 = sum((y - mean(y))^2)
if (flat) {
ish = ish + 1
		fit[ish, ] = 1:n*0 + mean(y)
		bhat[[ish]] = 0
		sse1 = sum((y - fit[ish, ])^2)
		if (infocrit == "BIC") {ic[s, ish] <- n * log(sse1) + log(n)}
		if (infocrit == "CIC") {ic[s, ish] <- log(sse1) + log(2 / (n - 1) + 1)}
}
## decreasing
if (dec) {
ish = ish + 1
		z = pmult0 %*% y
		ans = coneA(z, amat1, msg = msg)		
		fit[ish, ] = delta %*% uinv0 %*% ans$thetahat
		bhat[[ish]] = uinv0 %*% ans$thetahat
		sse2 = sum((y - fit[ish, ])^2)
		if (infocrit == "BIC") {ic[s, ish] <- n * log(sse2) + log(n) * edf0[ish]}
		if (infocrit == "CIC") {ic[s, ish] <- log(sse2) + log(2 * (edf0[ish] + 1) / (n - 1 - 1.5 * edf0[ish]) + 1)}
#new:
		if ((2 * (edf0[ish] + 1) / (n - 1 - 1.5 * edf0[ish]) + 1) < 0) {stop ("The sample size n is too small to make a valid CIC for a decreasing shape!")}
}
if (jp) {
ish = ish + 1
## jump --need to loop through possible jumppoints
		mj = length(delta) / n + 2
		minsse = sum((y - mean(y))^2)
#rm_vec = NULL
		for (i in 1:(n - 1)) {
			smat = matrix(0, nrow = k + 3, ncol = mj)
			smat[1:k, 1:(mj - 2)] = -slopes
			smat[k + 1, mj - 1] = 1

			djump = matrix(0, nrow = n, ncol = mj)
			djump[ ,1:(mj - 2)] = delta

			djump[1:i, mj - 1] = (2 * i + 1) / 2 - n; djump[(i + 1):n, mj - 1] = (2 * i + 1) / 2
			djump[1:i, mj] = 0; djump[(i + 1):n, mj] = x[x > (x[i] + x[i + 1]) / 2] - (x[i] + x[i + 1] ) / 2
			m_i0 = mean(djump[, mj])
			djump[, mj] = djump[, mj] - mean(djump[, mj])

			kn = 1:(k + 3) < 0
			kn[1:k] = knots > (x[i] + x[i + 1])/2
			smat[, mj] = 0
			smat[kn, mj] = -1; 
			smat[k + 2, mj] = -1
			ansi = -sl((x[i] + x[i + 1]) / 2, knots, slopes)
			smat[k + 2, 1:m0] = ansi 
			smat[k + 3, 1:m0] = ansi 

			use = 1:(k + 3) > 0
			kp = min(knots[knots > (x[i] + x[i + 1]) / 2])
			if (sum(x > (x[i] + x[i + 1]) / 2 & x <= kp) == 0) {use[k + 2] = FALSE}
			kp = max(knots[knots < (x[i] + x[i + 1]) / 2])
			if (sum(x < (x[i] + x[i + 1]) / 2 & x >= kp) == 0) {use[k + 3] = FALSE}
			smat = smat[use, ]

			if (i == 1 | i == (n - 1)) {
				smat = smat[, -mj, drop = FALSE]
				djump = djump[, -mj, drop = FALSE] 
			}
	
			#umat = chol(t(djump) %*% djump)
			umat = chol(crossprod(djump))
			uinv = solve(umat)
			pmult = t(uinv) %*% t(djump)
			z = pmult %*% y
			amat = smat %*% uinv
			fiti = coneA(z, amat, msg = msg)
			phat = fiti$thetahat
			theta = djump %*% uinv %*% phat
			bi = uinv %*% phat 
			ssei = sum((y - theta)^2)
			if (ssei < minsse) {ijp = i; sse3 = ssei; thb = theta; df = fiti$df; minsse = ssei; b3 = bi; m_i = m_i0} 
		}
		fit[ish, ] = thb
		bhat[[ish]] = b3
		ijps0[ish] = ijp
		m_is0[ish] = m_i
		if (infocrit == "BIC") {ic[s, ish] <- n * log(sse3) + log(n) * edf0[ish]}
		if (infocrit == "CIC") {ic[s, ish] <- log(sse3) + log(2 * (edf0[ish] + 1) / (n - 1 - 1.5 * edf0[ish]) + 1)}
#new:
		if ((2 * (edf0[ish] + 1) / (n - 1 - 1.5 * edf0[ish]) + 1) < 0) {stop ("The sample size n is too small to make a valid CIC for a one-jump shape!")}	
}	
if (invee) {
ish = ish + 1
## inverted vee -- need to loop through possible inter-knot spaces!!
		minsse = sum((y - mean(y))^2)
		av = slopes
		#cv = t(delta) %*% y
		cv = crossprod(delta, y)
		for (j in 2:k) {
		#for (j in 1:k) {		
			av1 = av
			av1[j:k,] = -av1[j:k,]
			qans = qprog(qv, cv, av1, 1:k*0, msg = msg)
			theta = delta %*% qans$thetahat
			bj = qans$thetahat 
			sse = sum((y - theta)^2)
			if (sse < minsse) {ch = j; thb = theta; minsse = sse; df = qans$df; sse4 = sse; b4 = bj; pos = k}
		}
		fit[ish, ] = thb
		bhat[[ish]] = b4
		ijps0[ish] = pos
#	edf4=findED2(delta,x,1000)
		if (infocrit == "BIC") {ic[s, ish] <- n * log(sse4) + log(n) * edf0[ish]}
		if (infocrit == "CIC") {ic[s, ish] <- log(sse4) + log(2 * (edf0[ish] + 1) / (n - 1 - 1.5 * edf0[ish]) + 1)}
#new:
		if ((2 * (edf0[ish] + 1) / (n - 1 - 1.5 * edf0[ish]) + 1) < 0) {stop ("The sample size n is too small to make a valid CIC for an inverted vee shape!")}
}
if (vee) {
ish = ish + 1
## vee -- need to loop through possible inter-knot spaces!!
		minsse = sum((y - mean(y))^2)
		av = slopes
		cv = crossprod(delta, y)
		for (j in 2:k) {
		#for (j in 1:k) {
			av1 = -av
			av1[j:k,] = -av1[j:k,]
			qans = qprog(qv, cv, av1, 1:k*0, msg = msg)
			theta = delta %*% qans$thetahat
			bj = qans$thetahat
			sse = sum((y - theta)^2)
			if (sse < minsse) {ch = j; thb = theta; minsse = sse; df = qans$df; sse5 = sse; b5 = bj; pos = k}
		}
		fit[ish, ] = thb
		bhat[[ish]] = b5
		ijps0[ish] = pos 
#	edf4=findED2(delta,x,1000)
		if (infocrit == "BIC") {ic[s, ish] <- n * log(sse5) + log(n) * edf0[ish]}
		if (infocrit == "CIC") {ic[s, ish] <- log(sse5) + log(2 * (edf0[ish] + 1) / (n - 1 - 1.5 * edf0[ish]) + 1)}
#new:
		if ((2 * (edf0[ish] + 1) / (n - 1 - 1.5 * edf0[ish]) + 1) < 0) {stop ("The sample size n is too small to make a valid CIC for a vee shape!")}
}
if (inc) {
ish = ish + 1
## increasing - ggm
		my.lm = lm((-y) ~ 1 + x)
		fit[ish,] = -my.lm$fitted.values
		#fit[6,] = -fit[6,]
		bhat[[ish]] = as.matrix(unname(coef(my.lm)))
	##consider only if positive slope - else just assign sse1+.01 - jerryrig
		if (my.lm$coefficients[2] < 0) {sse6 = sum((y - fit[ish, ])^2)}
		if (my.lm$coefficients[2] >= 0) {sse6 = sse1 + .01}

	#ggm
		if (infocrit == "BIC") {ic[s, ish] <- n * log(sse6) + log(n) * edf0[ish]}
		if (infocrit == "CIC") {ic[s, ish] <- log(sse6) + log(2 * (edf0[ish] + 1) / (n - 1 - 1.5 * edf0[ish]) + 1)}
#new:
		if ((2 * (edf0[ish] + 1) / (n - 1 - 1.5 * edf0[ish]) + 1) < 0) {stop ("The sample size n is too small to make a valid CIC for an increasing shape!")}
#new: db choice 
}
if (db) {
ish = ish + 1
# two jumps -- need to loop through possible jumppoints
		m0 = length(delta) / n
		mj =  m0 + 4
		minsse = sum((y - mean(y))^2)
		rm_vec = NULL
		for (i in 1:(n - 4)) {
			for (j in (i + 3):(n - 1)) {
				smat = matrix(0, nrow = k + 6, ncol = mj)
				smat[1:k, 1:(mj - 4)] = -slopes
				smat[k + 1, mj - 3] = 1
				smat[k + 4, mj - 1] = 1
#new:
				#smat[k + 4, mj - 3] = 1

				djump = matrix(0, nrow = n, ncol = mj)
				djump[ ,1:(mj - 4)] = delta

				djump[1:i, mj - 3] = (2 * i + 1) / 2 - n; djump[(i + 1):n, mj - 3] = (2 * i + 1) / 2
				djump[1:i, mj - 2] = 0; djump[(i + 1):n, mj - 2] = x[x > (x[i] + x[i + 1]) / 2] - (x[i] + x[i + 1] ) / 2
				m_i0 = mean(djump[, mj - 2])
				djump[, mj - 2] = djump[, mj - 2] - mean(djump[, mj - 2])
		
				djump[1:j, mj - 1] = (2 * j + 1) / 2 - n; djump[(j + 1):n, mj - 1] = (2 * j + 1) / 2
				djump[1:j, mj] = 0; djump[(j + 1):n, mj] = x[x > (x[j] + x[j + 1]) / 2] - (x[j] + x[j + 1]) / 2
				m_j0 = mean(djump[, mj])
				djump[, mj] = djump[, mj] - mean(djump[, mj])
			
				kn = 1:(k + 6) < 0
				kn[1:k] = knots > (x[i] + x[i + 1]) / 2			
				smat[kn, mj - 2] = -1 
				smat[k + 2, mj - 2] = -1
				ansi = -sl((x[i] + x[i+1]) / 2, knots, slopes)
				smat[k + 2, 1:m0] = ansi
				smat[k + 3, 1:m0] = ansi 

				kn = 1:(k + 6) < 0
				kn[1:k] = knots > (x[j] + x[j + 1]) / 2	
				smat[kn, mj] = -1

				smat[k + 5, mj] = -1
				smat[k + 5, mj - 2] = -1
				ansj = -sl((x[j] + x[j+1]) / 2, knots, slopes)
				smat[k + 5, 1:m0] = ansj 

				smat[k + 6, 1:m0] = ansj 
				smat[k + 6, mj - 2] = -1 

				use = 1:(k + 6) > 0
				kp = min(knots[knots > (x[i] + x[i + 1]) / 2])
				if (sum(x > (x[i] + x[i + 1]) / 2 & x <= kp) == 0){use[k + 2] = FALSE}
				kp = max(knots[knots < (x[i] + x[i + 1]) / 2 ])
				if (sum(x < (x[i] + x[i + 1]) / 2 & x >= kp) == 0){use[k + 3] = FALSE}
				kp = min(knots[knots > (x[j] + x[j + 1]) / 2])
				if (sum(x > (x[j] + x[j + 1]) / 2 & x <= kp) == 0){use[k + 5] = FALSE}
				kp = max(knots[knots < (x[j] + x[j + 1]) / 2])	
				if (sum(x < (x[j] + x[j + 1])/2 & x >= kp) == 0){use[k + 6] = FALSE}	
				smat = smat[use, ]
	
				if (i == 1) {
					rm_vec = c(rm_vec, mj - 2)
				}

				if (j == (n - 1)) {
					rm_vec = c(rm_vec, mj)
				}

				if (!is.null(rm_vec)) {
					smat = smat[, -rm_vec, drop = FALSE]
					djump = djump[, -rm_vec, drop = FALSE]
					rm_vec = NULL
				}
				#umat = chol(t(djump) %*% djump)
				umat = chol(crossprod(djump))
				uinv = solve(umat)
				pmult = t(uinv) %*% t(djump)
				z = pmult %*% y
				amat = smat %*% uinv			
				fiti = coneA(z, amat, msg = msg)
				phat = fiti$thetahat
				theta = djump %*% uinv %*% phat
				bi = uinv %*% phat
				ssei = sum((y - theta)^2)
				if (ssei < minsse) {ijp = i; jjp = j; sse7 = ssei; thb = theta; df = fiti$df; minsse = ssei; b7 = bi; m_i = m_i0; m_j = m_j0}
			}
		}
		fit[ish, ] = thb
		bhat[[ish]] = b7
		ijps0[ish] = ijp
		jjps0[ish] = jjp
		m_is0[ish] = m_i
		m_js0[ish] = m_j 
		if (infocrit == "BIC") {ic[s, ish] <- n * log(sse7) + log(n) * edf0[ish]}
#new:
		if ((2 * (edf0[ish] + 1) / (n - 1 - 1.5 * edf0[ish]) + 1) < 0) {stop ("The sample size n is too small to make a valid CIC for a double-jump shape!")}
		if (infocrit == "CIC") {ic[s, ish] <- log(sse7) + log(2 * (edf0[ish] + 1) / (n - 1 - 1.5 * edf0[ish]) + 1)}
}
### select best shape for pixel ###	
		sel = which(ic[s, ] == min(ic[s, ]))
		bshp = shape[sel]
		thetab[s, ] = fit[sel, ]
		shb[s] = bshp 
		bs[[s]] = bhat[[sel]]
		ijps[[s]] = ijps0[[sel]]
		jjps[[s]] = jjps0[[sel]]
		m_is[[s]] = m_is0[[sel]]
		m_js[[s]] = m_js0[[sel]]
	}
	ans = list(shape = shb, ic = t(ic), thetab = t(thetab), fit = t(fit), x = x0, ymat = ymat, infocrit = infocrit, k = k, bs = bs, ijps = ijps, jjps = jjps, m_is = m_is, m_js = m_js, shp_in = bool)
	class(ans) = "ShapeSelectForest"
	return (ans)
}



########################################
#       MAKE THE EDGE VECTORS          #
########################################
#######################
#bqspline basis matrix# 
#######################
bqspl = function(x, m, knots = NULL, pic = FALSE, spl = TRUE) {
############
#make knots#  
############
	bmat = slopes = NULL
	xu = unique(x)
#knots are equal x quantile of length = m spaced on the support of min(x) to max(x)
	if (is.null(knots)) {
		#knots = 0:(m - 1) / (m - 1) * (max(x) - min(x)) + min(x)
		knots = quantile(xu, probs = seq(0, 1, length = m), names = FALSE)
		t = 1:(m+4)*0
		t[1] = t[2] = min(x)
		t[m+3] = t[m+4] = max(x)
		t[3:(m+2)] = knots
	}
##############################
#make quadratic bspline basis#
##############################
	n = length(x)
if (spl) {
	bmat = matrix(0, nrow = n, ncol = (m+1))
#1st edge
	bool = x <= t[4] & x >= t[3]
	bmat[bool, 1] = (t[4] - x[bool]) / (t[4] - t[2]) * (t[4] - x[bool]) / (t[4] - t[3])

#2nd edge 1st part 
	bool = x <= t[4] & x >= t[3]
	bmat[bool, 2] = (x[bool] - t[2]) / (t[4] - t[2]) * (t[4] - x[bool]) / (t[4] - t[3]) + (t[5] - x[bool]) / (t[5] - t[3]) * (x[bool] - t[3]) / (t[4] - t[3]) 

#2nd edge 2nd part 
	bool = x <= t[5] & x >= t[4]
	bmat[bool, 2] = (t[5] - x[bool]) / (t[5] - t[3]) * (t[5] - x[bool]) / (t[5] - t[4])

#3rd edge to the (m-1)th edge
	for(i in 3:(m-1)) {
		#1st part
		bool = x <= t[i+1] & x >= t[i]
		bmat[bool, i] = (x[bool] - t[i])**2 / (t[i+2] - t[i]) / (t[i+1] - t[i])
	
		#2nd part
		bool = x <= t[i+2] & x >= t[i+1]
		bmat[bool, i] = (x[bool] - t[i]) * (t[i+2] - x[bool]) / (t[i+2] - t[i]) / (t[i+2] - t[i+1]) + (t[i+3] - x[bool]) * (x[bool] - t[i+1]) / (t[i+3] - t[i+1]) / (t[i+2] - t[i+1])

		#3rd part
		bool = x <= t[i+3] & x >= t[i+2]
		bmat[bool, i] = (t[i+3] - x[bool])**2 / (t[i+3] - t[i+1]) / (t[i+3] - t[i+2])
	}

#the mth edge 1st part
	bool = x <= t[m+1] & x >= t[m]
	bmat[bool, m] = (x[bool] - t[m])**2 / (t[m+2] - t[m]) / (t[m+1] - t[m])
	 
#the mth edge 2nd part 
	bool = x <= t[m+2] & x >= t[m+1]
	bmat[bool, m] = (x[bool] - t[m]) * (t[m+2] - x[bool]) / (t[m+2] - t[m]) / (t[m+2] - t[m+1]) + (t[m+3] - x[bool]) * (x[bool] - t[m+1]) / (t[m+3] - t[m+1]) / (t[m+2] - t[m+1])

#the (m+1)th edge 
	bool = x <= t[m+2] & x >= t[m+1]
	bmat[bool, m+1] = (x[bool] - t[m+1])**2 / (t[m+3] - t[m+1]) / (t[m+2] - t[m+1])

#####################################################################
#make bqsplines' 1st derivatives at knots. row: knots; column: basis#
#####################################################################
	slopes = matrix(0, nrow = m, ncol = (m+1))
#edge 1
	slopes[1,1] = -2 / (t[4] - t[2])
	#slopes[2,1] = 0 	
	
#edge 2 
	slopes[1,2] = 2 / (t[4] - t[2])
	slopes[2,2] = -2 / (t[5] - t[3])
	#slopes[3,2] = 0

#edge 3 ~ m-1
	for(i in 3:(m-1)){
		#slopes[i-2,i] = 0
		slopes[i-1,i] = 2 / (t[i+2] - t[i])
		slopes[i,i] = -2 / (t[i+3] - t[i+1])
		#slopes[i+1,i] = 0
	}

#edge m
	#slopes[m-2,m] = 0
	slopes[m-1,m] = 2 / (t[m+2] - t[m])
	slopes[m,m] = -2 / (t[m+3] - t[m+1])

#edge m+1
	#slopes[m-1,m+1] = 0
	slopes[m,m+1] = 2 / (t[m+3] - t[m+1])
}
#####################
#plot bqspline basis#
#####################
	xpl = bpl = NULL
if (pic) {
	xpl = 0:1000/1000*(max(x)-min(x))+min(x)
	bpl = matrix(1:(1001*(m+1))*0,nrow=1001)
	
#1st edge
	bool = xpl <= t[4] & xpl >= t[3]
	bpl[bool, 1] = (t[4] - xpl[bool]) / (t[4] - t[2]) * (t[4] - xpl[bool]) / (t[4] - t[3])

#2nd edge 1st part 
	bool = xpl <= t[4] & xpl >= t[3]
	bpl[bool, 2] = (xpl[bool] - t[2]) / (t[4] - t[2]) * (t[4] - xpl[bool]) / (t[4] - t[3]) + (t[5] - xpl[bool]) / (t[5] - t[3]) * (xpl[bool] - t[3]) / (t[4] - t[3]) 

#2nd edge 2nd part 
	bool = xpl <= t[5] & xpl >= t[4]
	bpl[bool, 2] = (t[5] - xpl[bool]) / (t[5] - t[3]) * (t[5] - xpl[bool]) / (t[5] - t[4])

#3rd edge to the (m-1)th edge
	for(i in 3:(m-1)){
		#1st part
		bool = xpl <= t[i+1] & xpl >= t[i]
		bpl[bool, i] = (xpl[bool] - t[i])**2 / (t[i+2] - t[i]) / (t[i+1] - t[i])
	
		#2nd part
		bool = xpl <= t[i+2] & xpl >= t[i+1]
		bpl[bool, i] = (xpl[bool] - t[i]) * (t[i+2] - xpl[bool]) / (t[i+2] - t[i]) / (t[i+2] - t[i+1]) + (t[i+3] - xpl[bool]) * (xpl[bool] - t[i+1]) / (t[i+3] - t[i+1]) / (t[i+2] - t[i+1])

		#3rd part
		bool = xpl <= t[i+3] & xpl >= t[i+2]
		bpl[bool, i] = (t[i+3] - xpl[bool])**2 / (t[i+3] - t[i+1]) / (t[i+3] - t[i+2])
	}

#the mth edge 1st part
	bool = xpl <= t[m+1] & xpl >= t[m]
	bpl[bool, m] = (xpl[bool] - t[m])**2 / (t[m+2] - t[m]) / (t[m+1] - t[m])
	 
#the mth edge 2nd part 

	bool = xpl <= t[m+2] & xpl >= t[m+1]
	bpl[bool, m] = (xpl[bool] - t[m]) * (t[m+2] - xpl[bool]) / (t[m+2] - t[m]) / (t[m+2] - t[m+1]) + (t[m+3] - xpl[bool]) * (xpl[bool] - t[m+1]) / (t[m+3] - t[m+1]) / (t[m+2] - t[m+1])

#the (m+1)th edge 
	bool = xpl <= t[m+2] & xpl >= t[m+1]
	bpl[bool, m+1] = (xpl[bool] - t[m+1])**2 / (t[m+3] - t[m+1]) / (t[m+2] - t[m+1])
}
	ans = new.env()
	ans$bmat = bmat
	ans$bpl = bpl
	ans$xpl = xpl 
	ans$slopes = slopes
	ans$knots = knots
	ans
}

############################################################
#make the first derivative at a jump point by interpolation#
############################################################
sl = function(xinterp, knots, slopes) {
	nk = length(knots)
	obs = 1:nk
	dist = knots - xinterp
	id = NULL; id1 = NULL
	if (any(round(dist, 8) == 0)) {
		id = min(obs[round(dist, 8) == 0])
	} else {	
		for (i in 1:(nk - 1)) {
			if (dist[i] < 0 & dist[i + 1] > 0) {
				id = i
				break 
			}
		}
	}
	k1 = knots[id]
	k2 = knots[id + 1]
	sinterp = 1:(nk + 1)*0
	a = (xinterp - k1) / (k2 - k1) 
	yk1 = slopes[id, ]
	yk2 = slopes[id + 1, ]
	sinterp = a * (yk2 - yk1) + yk1
	sinterp
}


#####################
#make a plot of fits#
##################### 
plotshape = function(object, ids = 1, color = "mediumorchid4", lty = 1, lwd = 1, cex = .83, cex.main = .93, form = TRUE, icpic = FALSE, both = TRUE, tt = NULL, transpose = FALSE, plot = graphics:::plot) {
	if (!inherits(object, "ShapeSelectForest")) { 
	        warning("calling plotpersp(<fake-ShapeSelectForest-object>) ...")
        }
	fits = object$thetab
	x0 = object$x; ymat = object$ymat; k = object$k; bs = object$bs; ic = object$ic; infocrit = object$infocrit
	shape = object$shape; shp_in = object$shp_in; ijps = object$ijps; jjps = object$jjps; m_is = object$m_is; m_js = object$m_js; 
	x0s = sort(x0)
	x = (x0s - min(x0s)) / (max(x0s) - min(x0s))
	n = ncol(ymat)
	nr = nrow(ymat)
	#shps = c('flat', 'decr', 'one-jp', 'invee', 'vee', 'incr', 'two-jp')
	shps = c('flat', 'decr', 'one-jp', 'invee', 'vee', 'incr', 'two-jp')[shp_in]
	xlst = ylst = flst = ynmlst = list()
	tt0 = tt
	ans0 = 	bqspl(x, k, pic = TRUE, spl = FALSE)
	xpl = ans0$xpl
	xpl0 = 0:1000/1000*(max(x0s) - min(x0s)) + min(x0s)
	bpl = ans0$bpl  
	nthps = length(bs)
	thps = list()
	for (i in 1:nthps) {
		b = bs[[i]]; sh = shape[i]
		if (sh == 1) {
			thps[[i]] = 1:1001*0 + mean(ymat[, i])
		} else if (sh == 6) {
			thps[[i]] = -cbind(1:1001*0 + 1, xpl) %*% b
		} else {
			if (sh == 3 | sh == 7) {
				ijp = ijps[[i]]; jjp = jjps[[i]]
				m_i = m_is[[i]]; m_j = m_js[[i]]
			
				d1 = d2 = 1:1001*0
				d1[xpl <= (x[ijp] + x[ijp + 1]) / 2] = (2 * ijp + 1) / 2 - nr
				d1[xpl > (x[ijp] + x[ijp + 1]) / 2] = (2 * ijp + 1) / 2
				d2[xpl <= (x[ijp] + x[ijp + 1]) / 2] = 0
				d2[xpl > (x[ijp] + x[ijp + 1]) / 2] = xpl[xpl > (x[ijp] + x[ijp+1]) / 2] - (x[ijp] + x[ijp+1]) / 2
				d2 = d2 - m_i
				#di = cbind(d1, d2)
				if (sh == 7) {
					d3 = d4 = 1:1001*0
					d3[xpl <= (x[jjp] + x[jjp+1]) / 2] = (2 * jjp + 1) / 2 - nr
					d3[xpl > (x[jjp] + x[jjp+1]) / 2] = (2 * jjp + 1) / 2
					d4[xpl <= (x[jjp] + x[jjp+1]) / 2] = 0
					d4[xpl > (x[jjp] + x[jjp+1]) / 2] = xpl[xpl > (x[jjp] + x[jjp+1]) / 2] - (x[jjp] + x[jjp+1]) / 2
					d4 = d4 - m_j 
					#di = cbind(di, d3, d4)
				}
				if (sh == 3) { 
					if (ijp == 1 | ijp == (nr - 1)) {d2 = NULL}
					di = cbind(d1, d2)
				}
				if (sh == 7) { 
					if (ijp == 1) {d2 = NULL}
					if (jjp == (nr - 1)) {d4 = NULL}
					di = cbind(d1, d2, d3, d4)
				}
				bpl1 = cbind(bpl, di)
			}  else {
				bpl1 = bpl
			}	
			thps[[i]] = bpl1 %*% b
		}
	}
	if (both) {
		k = 2
		xlst[[1]] = x0s; xlst[[2]] = 1:nrow(ic)
		ylst[[1]] = ymat; ylst[[2]] = ic
		flst[[1]] = fits; flst[[2]] = ic
		ynmlst[[1]] = 'y'; ynmlst[[2]] = 'ic'
	} else {
		k = 1
		if (icpic) {
			xlst[[1]] = 1:nrow(ic); ylst[[1]] = ic; flst[[1]] = ic; ynmlst[[1]] = 'ic'
			#x = 1:nrow(ic); ys = ic; fs = ic; ynm = 'ic'
		} else {
			xlst[[1]] = x0s; ylst[[1]] = ymat; flst[[1]] = fits; ynmlst[[1]] = 'y'
			#ys = ymat; fs = fits; ynm = 'y'
		}
	}
	palette = colors()[c(552, 498, 654, 254, 635, 26, 32, 547)] 
	if (!all(ids %in% (1:n))) {
		stop ('All ids must be within the range of 1 to the number of columns of ymat!')
	} 
	nids0 = length(ids)
#new: turn off both	
if (icpic) {k = 1; xlst[[1]] = 1:nrow(ic); ylst[[1]] = ic; flst[[1]] = ic; ynmlst[[1]] = 'ic'}
	nids = nids0 * k
	if (form) {
		width = nids^.5
		wd = round(width)
		if (width > wd) {
			wd = wd + 1
		}  
		if ((wd^2 - nids) >= wd ) {
			fm = c(wd, wd - 1)
		} else {
			fm = rep(wd, 2)
		}
		if (wd > 3) {
			cex = .7
			cex.main = .8
		}
		if (transpose) {
			fm = rev(fm)
		}
	}
	#if (static) {
	if (form) {
		par(mfrow = c(fm[1], fm[2]))
		par(mar = c(4, 1, 1, 1))
		par(cex.main = cex.main)
	}
	for (i in 1:nids0) {		
		ch = ids[i] 	
		for (j in 1:k) {
			xs = xlst[[j]]; ys = ylst[[j]]; fs = flst[[j]]; ynm = ynmlst[[j]]
			if (j == 1) {
				if (!icpic) {
					r1 = range(ys[ ,ch]); r2 = range(thps[[ch]]); r3 = range(fs[ ,ch])
					rmin = min(c(r1, r2, r3)); rmax = max(c(r1, r2, r3))
					rg = rmax - rmin 
					plot(xs, ys[ ,ch], col = 1, cex = 1.1, ylab = ynm, xlab = '', ylim = c(rmin, rmax))
					lines(xpl0, thps[[ch]], type = 'l', lty = 2, xlab = '', ylab = ynm, col = color, lwd = 2)
					points(xs, fs[ ,ch], col = color, pch = 20, cex = cex) 
					#legend('topleft', bty = 'n', expression(hat(f)), col = color, lty = 2, lwd = 2)
				} else {
					plot(xs, ys[ ,ch], type = 'l', xaxt = 'n', ann = TRUE, xlab = '', ylab = ynm)
					points(xs, fs[ ,ch], col = color, pch = 20, cex = 1.3)
					abline(h = min(ys[ ,ch]), lty = 2, lwd = 2) 
					ticks = shps[1:nrow(ic)]
					bool = 1:nrow(ic) %% 2 != 0
					mtext(ticks[bool], side = 1, at = (1:nrow(ic))[bool], line = 1, cex = .83)	
					mtext(ticks[!bool], side = 1, at = (1:nrow(ic))[!bool], line = .1,  cex = .83)	
				}
			} else if (j == 2) {
				plot(xs, ys[ ,ch], type = 'l', xaxt = 'n', ann = TRUE, xlab = '', ylab = ynm)
				points(xs, fs[ ,ch], col = color, pch = 20, cex = 1.3)
				abline(h = min(ys[ ,ch]), lty = 2, lwd = 2) 
				ticks = shps[1:nrow(ic)]
				if (nrow(ic) > 2) {
					bool = 1:nrow(ic) %% 2 != 0
					mtext(ticks[bool], side = 1, at = (1:nrow(ic))[bool], line = 1, cex = .83)	
					mtext(ticks[!bool], side = 1, at = (1:nrow(ic))[!bool], line = .1,  cex = .83)	
				} else if (nrow(ic) == 2) { 
					bool = 1:nrow(ic) %% 2 != 0
					mtext(ticks[bool], side = 1, at = (1:nrow(ic))[bool], line = .1, cex = .83)	
					mtext(ticks[!bool], side = 1, at = (1:nrow(ic))[!bool], line = .1,  cex = .83)
				} else {
					mtext(ticks, side = 1, at = (1:nrow(ic)), line = .1, cex = .83)
				}
			}
			if (is.null(tt0)) {	
				if (j == 2 | icpic) {
					if (ch == 1) {
						tt = paste(paste(infocrit, "'s for"), 'the 1st Col of Ymat') 
					} else if (ch == 2) {
						tt = paste(paste(infocrit,  "'s for"), 'the 2nd Col of Ymat') 
					} else if (ch == 3) {
						tt = paste(paste(infocrit,  "'s for"), 'the 3rd Col of Ymat') 
					} else {
						tt = paste(paste(infocrit,  "'s for"), paste('the', paste(ch,'th Col of Ymat', sep = ''))) 
					}
				} else {
					if (ch == 1) {
						#tt = 'The 1st col of ymat'
						tt = paste('The Best Shape for', 'the 1st Col of Ymat') 
					} else if (ch == 2) {
						#tt = 'The 2nd col of ymat'
						tt = paste('The Best Shape for', 'the 2nd Col of Ymat') 
					} else if (ch == 3) {
						#tt = 'The 3rd col of ymat'
						tt = paste('The Best Shape for', 'the 3rd Col of Ymat') 
					} else {
						#tt = paste('The', paste(ch,'th col of ymat', sep = ''))
						tt = paste('The Best Shape for', paste('the', paste(ch, 'th Col of Ymat', sep = ''))) 
					}
				}
			} else {
				if (j == 2) {
					tt = paste(paste(infocrit, "'s for"), tt0[i]) 
				} else {
					tt = paste('The Best Shape for', tt0[i]) 
				}					
			}
			title(tt)
			#lines(x, fs[ ,ch], col = color, lty = lty, lwd = lwd)
		}
	}
	#} #else {
	#	for (i in 1:nids0) {
	#		par(mfrow = c(1, 1))
	#		ch = ids[i] 
	#		for (j in 1:k) {
	#			x = xlst[[j]]; ys = ylst[[j]]; fs = flst[[j]]; ynm = ynmlst[[j]]
	#			if (j == 1 & !icpic) {
	#				plot(x, ys[ ,ch], ylab = ynm)
	#			} else if (j == 2 | icpic) {
	#				plot(x, ys[ ,ch], xaxt = 'n', ann = TRUE, ylab = ynm)
	#			}
	#			if (j == 2 | icpic) {
	#				abline(h = min(ys[ ,ch]), lty = 2) 
	#				ticks = shps[1:nrow(ic)]
	#				bool = 1:nrow(ic) %% 2 != 0
	#				mtext(ticks[bool], side = 1, at = (1:nrow(ic))[bool], line = 1, cex = cex)	
	#				mtext(ticks[!bool], side = 1, at = (1:nrow(ic))[!bool], line = .1, cex = cex)	
	#			}
	#			if (is.null(tt0)) {	
	#				if (j == 2 | icpic) {
	#					if (ch == 1) {
	#						tt = paste(paste(infocrit, "'s for"), 'the 1st Col of Ymat') 
	#					} else if (ch == 2) {
	#						tt = paste(paste(infocrit, "'s for"), 'the 2nd Col of Ymat') 
	#					} else if (ch == 3) {
	#						tt = paste(paste(infocrit, "'s for"), 'the 3rd Col of Ymat') 
	#					} else {
	#						tt = paste(paste(infocrit, "'s for"), paste('the', paste(ch, 'th Col of Ymat', sep = ''))) 
	#					}
	#				} else {
	#					if (ch == 1) {
	#						#tt = 'The 1st col of ymat'
	#						tt = paste('The Best Shape for', 'the 1st Col of Ymat') 
	#					} else if (ch == 2) {
	#						#tt = 'The 2nd col of ymat'
	#						tt = paste('The Best Shape for', 'the 2nd Col of Ymat') 
	#					} else if (ch == 3) {
	#						#tt = 'The 3rd col of ymat'
	#						tt = paste('The Best Shape for', 'the 3rd Col of Ymat') 
	#					} else {
	#						#tt = paste('The', paste(ch,'th col of ymat', sep = ''))
	#						tt = paste('The Best Shape for', paste('the', paste(ch,'th Col of Ymat', sep = ''))) 
	#					}
	#				}
	#			} else {
	#				if (j == 2 | icpic) {
	#					tt = paste(paste(infocrit, "'s for"), tt0[i]) 
	#				} else {
	#					tt = paste('The Best Shape for', tt0[i]) 
	#				}					
	#			}
	#			title(tt)
	#			lines(x, fs[ ,ch], col = color, lty = lty, lwd = lwd)
	#			if ((i * j) < nids) {
	#				Sys.sleep(interv)
	#			}	
	#		}
	#	}
	#}
}


###########################################################################################################
#### shapeparams : function to read results from shape.fun and extract parameters from the model fit.
####		Function set up to run on one pixels at a time.
####		Inputs include: shapenum,ic,and thetab output from the shape function as well as
####			a vetor of yrs
####		Function assumes there is one predicted value per year
############################################################################################################
shapeparams = function(shapenum, ic, thetab, x) {
	n = length(x)
	cnt = 1:n
	start.yr = x[1]
	end.yr = x[n]
	flag = 0
	my.min = min(thetab)
	my.max = max(thetab)  
	
	pre.rate = 0
	dist.yr = 0
	dist.mag = 0
	dist.mag2 = 0
	dist.dur = 0
	post.rate = 0
	
	pre2.rate = 0
	dist2.yr = 0	
	dist2.mag = 0
	dist2.mag2 = 0
	dist2.dur = 0
	post2.rate = 0

	# flat
	if (shapenum == 1) {
		shape = "flat"
	}

  	# decrease
	if (shapenum == 2) {
		shape = "decr"
		pre.rate = ((thetab[1] - thetab[n]) / n) / abs(thetab[1])
		post.rate = pre.rate
	}

  	# jump
	if (shapenum == 3) {
		shape = "jump"
		jumpstart = round(thetab - c(thetab[-1], thetab[n]), 2)
	  # if jumps are drops at beginning of series...call decreasing
		if (sum(jumpstart < 0) == 0) {
			temp = 1
			shapenum = 2
			shape = "decr"
			pre.rate = ((thetab[1] - thetab[n]) / n) / abs(thetab[1])
			post.rate = pre.rate
			flag = 1
		}
	  # for normal jumps
		if (sum(jumpstart < 0) != 0) {
			temp = 2
			jumpstartpos = min(cnt[jumpstart < 0])
			jumpendpos = min(cnt[(jumpstart >= 0) & (cnt > jumpstartpos)])
			dist.yr = x[jumpstartpos + 1]
			if (jumpstartpos == 1) {pre.rate = 0}
			if (jumpstartpos > 1) {
				pre.rate = (thetab[1] - thetab[jumpstartpos]) / (((dist.yr - 1) - start.yr) * abs(thetab[1]))
			}
			dist.mag = thetab[jumpendpos] - thetab[jumpstartpos]
			dist.mag2 = dist.mag / my.max
			dist.dur = jumpendpos - jumpstartpos
			post.rate = (thetab[jumpendpos] - thetab[n]) / ((end.yr - (dist.yr-1)) * abs(thetab[jumpendpos]))
		}
		
	}
  	# inv 
	if (shapenum == 4) {
		shape = "inv"
		jumpstart = round(thetab - c(thetab[-1], thetab[n]), 2)
		## bandaid, if the inv never decreases, then set an incr shape
		if (sum(jumpstart > 0) == 0) {
			shapenum = 6
			flag = 2 
		}
		## else continue with inv
		if (sum(jumpstart > 0) > 0) {	
			jumpstartpos = min(cnt[jumpstart < 0])
			jumpendpos = min(cnt[(jumpstart > 0) & (cnt > jumpstartpos)])
			temp.mag = thetab[jumpendpos] - thetab[jumpstartpos]
			temp.dur = jumpendpos - jumpstartpos
	## primary parameters from point of acceleration till recovery		
			stdrate = -temp.mag / temp.dur
			accelpos = min(cnt[jumpstart < stdrate])
		# bandaid
			if ((accelpos < 1) || (accelpos > jumpendpos)) {
				accelpos = jumpstartpos
			}
			if (is.na(accelpos)) {
				accelpos = jumpstartpos
			}
			dist.yr = x[accelpos]
			dist.dur = jumpendpos - accelpos
			pre.rate = 0
			post.rate = (thetab[jumpendpos] - thetab[n]) / ((end.yr - (dist.yr - 1)) * abs(thetab[jumpendpos]))
			dist.mag = thetab[jumpendpos] - thetab[accelpos]
			dist.mag2 = dist.mag / my.max
			dist.mag = dist.mag

	## secondary parameters from point of inflection till point 
	##	of acceleration
			if (jumpstartpos < accelpos) {		
				dist2.yr = x[jumpstartpos]
				dist2.dur = accelpos - jumpstartpos
				dist2.mag = thetab[accelpos] - thetab[jumpstartpos]
				dist2.mag2 = (dist2.mag / my.max)
				dist2.mag = dist2.mag
			}
		## end continue
		}	
	}
  	# vee... if inflecpos = 1 call increasing and leave nonparametric fit
	if (shapenum == 5) {
		shape = "vee"
		jumpstart = round(thetab - c(thetab[-1], thetab[n]), 2)
		jumpstartpos = min(cnt[jumpstart < 0])
		if (jumpstartpos == 1) {
			shapenum = 6
			flag = 2
		}
		temp.mag = thetab[n] - thetab[jumpstartpos]
		temp.dur = n - jumpstartpos
		if (jumpstartpos > 1) {
		## primary parameters from point of acceleration till recovery
		## except prerate which is prior to point of inflection		
			stdrate = -temp.mag / temp.dur
			accelpos = min(cnt[jumpstart < stdrate])
			if (is.na(accelpos)) {
				accelpos = jumpstartpos
			}
			dist.yr = x[accelpos]
			dist.dur = n - accelpos
			pre.rate = (thetab[1] - thetab[jumpstartpos]) / (((x[jumpstartpos] - 1) - start.yr) * abs(thetab[1]))
			post.rate = 0
			dist.mag = thetab[n] - thetab[accelpos]
			dist.mag2 = dist.mag / my.max
			dist.mag = dist.mag

	## secondary parameters from point of inflection till point 
	##	of acceleration
			if (jumpstartpos < accelpos) {		
				dist2.yr = x[jumpstartpos]
				dist2.dur = accelpos - jumpstartpos
				dist2.mag = thetab[accelpos] - thetab[jumpstartpos]
				dist2.mag2 = dist2.mag / my.max
				dist2.mag = dist2.mag
			}	
		}
	}

  	# increase
	if (shapenum == 6) {
		shape = "incr"
		pre.rate = 0
		dist.yr = start.yr
		dist.mag = thetab[n] - thetab[1]
		dist.mag2 = dist.mag / my.max
		dist.mag = dist.mag
		dist.dur = n
		flag = max(flag, 0, na.rm = TRUE)
	}

	# 2jump: note that primary parameters reflect highest magnitude jump,
	#		seconday parameters reflect lesser magnitude jump. Also note
	#		pre rates are the same, that is prior to the first jump since
	#		post rate on first jump is really pre of second

	if (shapenum == 7) {
		shape = "2jump"
		jumpstart = round(thetab - c(thetab[-1], thetab[n]), 2)
		jumpstartpos = min(cnt[jumpstart < 0])
		jumpendpos = min(cnt[(jumpstart >= 0) & (cnt > jumpstartpos)])
		jumpstartpos2 = max(cnt[jumpstart < 0])
		jumpendpos2 = min(cnt[(jumpstart >= 0) & (cnt > jumpstartpos2)])
		dista.yr = x[jumpstartpos + 1]
		distb.yr = x[jumpstartpos2 + 1]

		## rest of parameters for first jump
		prea.rate = (thetab[1] - thetab[jumpstartpos]) / (((dista.yr - 1) - start.yr) * abs(thetab[1]))
		dista.mag = thetab[jumpendpos] - thetab[jumpstartpos]
		dista.mag2 = abs(dista.mag / my.max)
		dista.dur = jumpendpos - jumpstartpos
		posta.rate = (thetab[jumpendpos] - thetab[jumpstartpos2]) / (((distb.yr - 1) - (dista.yr - 1)) * abs(thetab[jumpendpos]))
	
		## rest of parameters for second jump
		preb.rate = (thetab[1] - thetab[jumpstartpos2]) / (((dista.yr - 1) - start.yr) * abs(thetab[1]))
		distb.mag = thetab[jumpendpos2] - thetab[jumpstartpos2]
		distb.mag2 = abs(distb.mag / my.max)
		distb.dur = jumpendpos2 - jumpstartpos2
		postb.rate = (thetab[jumpendpos2] - thetab[n]) / ((end.yr - (distb.yr - 1)) * abs(thetab[jumpendpos2]))

		## assign primary and secondary parameters based on mag
		if (dista.mag >= distb.mag){
			dist.yr = dista.yr
			pre.rate = prea.rate
			dist.mag = dista.mag
			dist.mag2 = dista.mag2
			dist.dur = dista.dur
			post.rate = posta.rate
 
			dist2.yr = distb.yr
			pre2.rate = posta.rate
			dist2.mag = distb.mag
			dist2.mag2 = distb.mag2
			dist2.dur = distb.dur
			post2.rate = postb.rate
		} 
		if (dista.mag < distb.mag){
			dist.yr = distb.yr
			pre.rate = posta.rate
			dist.mag = distb.mag
			dist.mag2 = distb.mag2
			dist.dur = distb.dur
			post.rate = postb.rate
 
			dist2.yr = dista.yr
			pre2.rate = prea.rate
			dist2.mag = dista.mag
			dist2.mag2 = dista.mag2
			dist2.dur = dista.dur
			post2.rate = posta.rate
		} 
			flag = 0
	}

	my.ic = ic[shapenum]

	shapeparams = data.frame(shapenum, pre.rate, dist.yr, dist.mag, dist.mag2, dist.dur, post.rate, my.ic,flag, pre2.rate, dist2.yr, dist2.mag, dist2.mag2, dist2.dur, post2.rate)

## Mop up: 
	## missing values set to 0 and flag 6
	myfun = function(mydat) {sum(is.na(mydat)) > 0}
	crit = apply(shapeparams, 1, myfun)
	shapeparams$flag[crit] = 6
	shapeparams[is.na(shapeparams)] = 0

	## shapenum != to 1-7 sets shape to flat and zero for all other params
	##		and flag 3
	crit = is.na(match(shapeparams$shapenum, 1:7))
	shapeparams[crit, -1] = 0
	shapeparams$flag[crit] = 3 
	shapeparams$shapenum[crit] = 1

	## substantive negative values (<-.001) set to 0 and flag 4
	myfun = function(mydat) {sum(mydat < (-.001)) > 0}
	crit = apply(shapeparams, 1, myfun)
	shapeparams$flag[crit] = 4
	shapeparams[shapeparams < (-.001)] = 0

	## non substantive negative values set to 0 and don't flag
	shapeparams[shapeparams < 0] = 0

	## infinity values set to 0 and flag 5
	myfun = function(mydat) {sum(mydat == Inf) > 0}
	crit = apply(shapeparams, 1, myfun)
	shapeparams$flag[crit] = 5
	shapeparams[shapeparams == Inf] = 0

	return(shapeparams)
}




########################## internal functions ########################################

# f2a.median
# f2a.mean
# f2a.wrapper
# grd2gri

##################### custom version of median and mean ##############################

f2a.median <- function(x) {
    switch(length(x), x, mean(x[1]), sort(x)[2], mean(sort(x)[2:3]), sort(x)[3])
}

f2a.mean <- function(x) {
	sum(x) / length(x)
}

############################ subfunctions for rasters ###############################

f2a.wrapper <- function(x, years, N.bands, mtbs) {
	f2a.out <- flat2annual(years = years, flat.pred = x[1], all.shapes = x[(1:N.bands) + 1], all.dyrs = x[(1:N.bands) + N.bands + 1], all.durs = x[(1:N.bands) + 2 * N.bands + 1], mtbs = mtbs)
	return(f2a.out)
}


f2p.wrapper <- function(x, years, N.bands, mtbs) {
	f2p.out <- flat2parameter(years = years, flat.pred = x[1], all.shapes = x[(1:N.bands) + 1], all.dyrs = x[(1:N.bands) + N.bands + 1], 
all.durs = x[(1:N.bands) + 2 * N.bands + 1], all.mags = x[(1:N.bands) + 3 * N.bands + 1], mtbs = mtbs)
	return(f2p.out)
}


grd2gri <- function(x) {
	paste(strsplit(x, ".grd")[[1]], ".gri", sep = "")
}
			
f2a.map.jpeg <- function(years, folder, OUTPUT.fn, height = 10, width = 10 * (dim(mapgrid.dist)[2] / dim(mapgrid.dist)[1]), units = "in", res = 400) {
### Define Color Ramp ###
	MAP.CODES <- data.frame(row = 1:7, integercode = 0:6)
	MAP.CODES$colors <- c("white", "grey50", "red", "blue", "green4", "brown", "springgreen")
	MAP.CODES$names <- c("Unclassified", "Conversion", "Fire", "Harvest", "Stable", "Stress", "Recovery")
### Define years ###
	N.years <- length(years)
###Input image file of maps###
	mapgrid.dist <- brick(paste(folder, OUTPUT.fn, sep = "/"))
### output filename will have path, base and extension in every case
	OUTPUTbase  <- basename(OUTPUT.fn)				
#name and extension, no path
	OUTPUTsplit <- strsplit(OUTPUTbase, split = "\\.")[[1]]
	OUTPUTname <- OUTPUTsplit[1]  				
#name, no extension or path
	OUTPUText <- paste(".", OUTPUTsplit[2], sep = "")
#just extension
	OUTPUTpath <- dirname(OUTPUT.fn)				
#just path
	if (OUTPUTpath == ".") { 
		OUTPUTpath <- folder 
	}
	OUTPUTfn.noext <- file.path(OUTPUTpath, OUTPUTname)	
#path and name, no extension
	for (i in 1:N.years) {
		print(paste("year:",years[i]))
		integergrid <- mapgrid.dist[[i]]
		v <- getValues(mapgrid.dist[[i]])
		v <- MAP.CODES$row[match(v,MAP.CODES$integercode)]
		integergrid <- setValues(integergrid, v)
		jpeg(filename = paste(OUTPUTfn.noext, "-Disturbance-", years[i], ".jpg", sep = ""), height = height, width = width, 
		units = units, res = res)
#dev.new(width = 8, height = 10)
		opar <- par(mar = c(0, 0, 0, 0), xpd = NA)
		image(integergrid, col = MAP.CODES$colors, xlab = "", ylab = "", xaxt = "n", yaxt = "n", zlim = c(1, nrow(MAP.CODES)),
		main = "", asp = 1, bty = "n")
#mtext(paste(years[i],"                      "),side=1,line=-1.2,cex=3)
#mtext(years[i],side=3,line=-5,cex=3.5)
		legend("topright", inset = 0.03, legend = MAP.CODES$names, fill = MAP.CODES$colors, bg = "white", title = years[i], cex = 1.5)
		dev.off()
	}
}


######################### apply f2a to rasters ###########################

f2a.raster <- function(years, folder.in, folder.out, OUTPUT.fn, flat.pred.fn, INPUT.bands, layer.shape = 1, layer.dyr = 2, layer.dur = 5) {
	N.years <- length(years)
	N.bands <- length(INPUT.bands)
### Create filename for native raster format map output
	TMPfn.map <- rasterTmpFile(prefix = paste("raster_tmp_", OUTPUT.fn, "_map_", sep = ""))
### Create filename for predictor raster brick
	OUTPUTfn.brick <- rasterTmpFile(prefix = paste("raster_tmp_", OUTPUT.fn, "_brick_", sep = ""))
### Set data type
	data.type <- "INT1U"
	Ylev <- 0:7
### Build raster stack
	RAST <- vector("list", 0)
	RAST[[1]] <- raster(paste(folder.in, flat.pred.fn, sep = "/"), band = 1)
	for (b in 1:N.bands) {
		RAST[[b + 1]] <- raster(paste(folder.in, INPUT.bands[b], sep = "/"), band = layer.shape)
		RAST[[b + N.bands + 1]] <- raster(paste(folder.in, INPUT.bands[b], sep = "/"), band = layer.dyr)
		RAST[[b + 2 * N.bands + 1]] <- raster(paste(folder.in, INPUT.bands[b], sep = "/"), band = layer.dur)
	}
	RS <- stack(RAST)
	RB <- brick(RS, values = TRUE, filename = OUTPUTfn.brick, overwrite = TRUE)
	names(RB) <- c(flat.pred.fn, paste(INPUT.bands, "shape", sep = "_"), paste(INPUT.bands, "dyr", sep = "_"), paste(INPUT.bands, "dur", sep = "_"))
	print("brick done")
### Loop through rows
	out <- brick(RAST[[1]], nl = N.years, values = FALSE)
	names(out) <- paste("Y", years, sep = "")
	dataType(out) <- data.type
#NAvalue(out) <- -9999
	NAvalue(out) <- 0
	print("write start")
	print(OUTPUT.fn)
	print(paste("datatype =", data.type))
	out <- writeStart(out, filename = TMPfn.map, overwrite = TRUE, datatype = data.type)
	print("starting loops")
	for (r in 1:(dim(RB)[1])) {
		print(paste("rows =", r))
		v <- data.frame(getValues(RB, r))
	### deal with -9999 ###
		nonPredict <- apply(((v == -9999) | is.na(v)), 1, any)
		v.pred <- matrix(NA, nrow = nrow(v), ncol = N.years)
		colnames(v.pred) <- paste("Y", years, sep = "")
		if (any(!nonPredict)) {
			v.pred[!nonPredict, ] <- t(apply(v[!nonPredict, ], 1, f2a.wrapper, years = years, N.bands = N.bands))}
			writeValues(out, v.pred, r)	
		}
		out <- writeStop(out)
	#out<-setMinMax(out)
		writeRaster(out, paste(folder.out, OUTPUT.fn, sep = "/"), overwrite = TRUE, datatype = data.type)
		file.remove(TMPfn.map)
		file.remove(grd2gri(TMPfn.map))
		file.remove(OUTPUTfn.brick)
		file.remove(grd2gri(OUTPUTfn.brick))
	}
######################### apply f2p to rasters ###########################

f2p.raster <- function(years, folder.in, folder.out, OUTPUT.fn, flat.pred.fn, INPUT.bands, layer.shape = 1, layer.dyr = 2, layer.dur = 5, layer.mag = 3) {
	N.bands <- length(INPUT.bands)
### Create filename for native raster format map output
	TMPfn.map <- rasterTmpFile(prefix = paste("raster_tmp_", OUTPUT.fn, "_map_", sep = ""))
### Create filename for predictor raster brick
	OUTPUTfn.brick <- rasterTmpFile(prefix = paste("raster_tmp_", OUTPUT.fn, "_brick_" , sep = ""))
### Set data type
	data.type <- "FLT4S"
### Build raster stack
	RAST <- vector("list", 0)
	RAST[[1]] <- raster(paste(folder.in, flat.pred.fn, sep = "/"), band = 1)
	for (b in 1:N.bands) {
		RAST[[b+1]]  <- raster(paste(folder.in, INPUT.bands[b], sep = "/"), band = layer.shape)
	RAST[[b+N.bands+1]]  <- raster(paste(folder.in, INPUT.bands[b], sep = "/"), band = layer.dyr)
	RAST[[b+2*N.bands+1]] <- raster(paste(folder.in, INPUT.bands[b], sep = "/"), band = layer.dur)
		RAST[[b+3*N.bands+1]] <- raster(paste(folder.in, INPUT.bands[b], sep = "/"), band = layer.mag)
	}
	RS <- stack(RAST)
	RB <- brick(RS, values = TRUE, filename = OUTPUTfn.brick, overwrite = TRUE)
	names(RB) <- c(flat.pred.fn, paste(INPUT.bands, "shape", sep = "_"), paste(INPUT.bands, "dyr", sep = "_"), paste(INPUT.bands, "dur", sep = "_"), paste(INPUT.bands, "mag", sep = "_"))
	print("brick done")
### Loop through rows
	out <- brick(RAST[[1]], nl = (3 + N.bands), values = FALSE)
	names(out) <- c("flat.pred", "ann.yr", "ann.dur", paste("ann.mag.", 1:N.bands, sep = "") )
	dataType(out) <- data.type
	NAvalue(out) <- -9999
	print("write start")
	print(OUTPUT.fn)
	print(paste("datatype =", data.type))
	out <- writeStart(out, filename = TMPfn.map, overwrite = TRUE, datatype = data.type)
	print("starting loops")
	for (r in 1:(dim(RB)[1])) {
		print(paste("rows =", r))
		v <- data.frame(getValues(RB, r))
	### deal with -9999 ###
		nonPredict <- apply(((v == -9999)|is.na(v)), 1, any)
		v.pred <- matrix(NA, nrow = nrow(v), ncol = 3 + N.bands)
		colnames(v.pred) <- c("flat.pred", "ann.yr", "ann.dur", paste("ann.mag.", 1:N.bands, sep = "") )
		if (any(!nonPredict)) {
			v.pred[!nonPredict, ] <- t(apply(v[!nonPredict, ], 1, f2p.wrapper, years = years, N.bands = N.bands))
		}
		writeValues(out, v.pred, r)	
}

	out <- writeStop(out)
#out<-setMinMax(out)
	writeRaster(out, paste(folder.out, OUTPUT.fn, sep = "/"), overwrite = TRUE, datatype = data.type)
	file.remove(TMPfn.map)
	file.remove(grd2gri(TMPfn.map))
	file.remove(OUTPUTfn.brick)
	file.remove(grd2gri(OUTPUTfn.brick))
}

#############pixel based function - equivalent of flat2annual################################

flat2annual <- function(years, all.shapes, all.durs, all.dyrs, mtbs, flat.pred){
	n = length(years)
	cnt = 1:n
	#if not disturbed, set all years to flat.pred
	ann.yr = 0
	ann.dur = 0
	ann.pred = rep(flat.pred,n)  
	#if disturbed
	if ((flat.pred != 0) & (flat.pred != 4) & (flat.pred != 6) ) {
		ann.pred = rep(4, n)
		#define which shapes are used for median (ie all shapes not flat or decreasing)
		use.for.med <- (all.shapes >= 3)
		#if not stress and jumps are present, only use jumps for median
		if ((flat.pred != 5) & any(all.shapes == 3)) {
			use.for.med <- (all.shapes == 3)
		}
		N.med <- sum(use.for.med)
		if (N.med > 0) {
			if (N.med <= 5) {
				ann.yr = f2a.median(all.dyrs[use.for.med])
				ann.dur = f2a.median(all.durs[use.for.med])
			} else {
				ann.yr = median(all.dyrs[use.for.med])
				ann.dur = median(all.durs[use.for.med])}

			#ann.end = ann.yr + ann.dur - 1 #if I subtract 1 then all disturbances turn into recovery in 2010
				ann.end = ann.yr + ann.dur
				if (ann.end > max(years)) {ann.end <- max(years)}
			#rounds to the nearest year found in years
				STARTdiff <- years - ann.yr
				ENDdiff <- years - ann.end
				startpos = cnt[STARTdiff == min(STARTdiff[STARTdiff >= 0])] #round up	
				endpos = cnt[ENDdiff == min(ENDdiff[ENDdiff >= 0])] #round up
			#endpos   = cnt[ENDdiff   == max(ENDdiff  [ENDdiff  <=0])] #round down

			#if not conversion, set recovery values
				if (flat.pred != 1) {ann.pred[startpos:n] = 6}
				ann.pred[startpos:endpos] = flat.pred
		}
	}
	return(ann.pred)
}

#############pixel based function - equivalent of flat2annual################################

flat2parameter <- function(years, all.shapes, all.durs, all.dyrs, all.mags, mtbs, flat.pred) {
	n = length(years)
	cnt = 1:n
	n.bands <- length(all.shapes)
	#if not disturbed, set all years to flat.pred
	OUT <- c(flat.pred, rep(0, 2 + n.bands))
	#names(OUT)<-c("flat.pred","ann.yr","ann.dur","ann.mag.b5","ann.mag.fi","ann.mag.nbr","ann.mag.ndvi")
	#if disturbed
	if ((flat.pred != 0) & (flat.pred != 4) & (flat.pred != 6) ) {
		#define which shapes are used for median (ie all shapes not flat or decreasing)
		use.for.med <- (all.shapes >= 3)
		#if not stress and jumps are present, only use jumps for median
		if ((flat.pred != 5) & any(all.shapes == 3)) {
			use.for.med <- (all.shapes == 3)
		}
		N.med <- sum(use.for.med)
		if (N.med > 0) {
			if (N.med <= 5) {
				OUT[1] = flat.pred	
				Ydiff = years - f2a.median(all.dyrs[use.for.med])
				OUT[2] = years[Ydiff == min(Ydiff[Ydiff >= 0])]
				OUT[3] = f2a.median(all.durs[use.for.med])
				OUT[4:(3 + n.bands)] = all.mags	
			} else {
				OUT[1] = flat.pred	
				Ydiff = years - median(all.dyrs[use.for.med])
				OUT[2] = years[Ydiff == min(Ydiff[Ydiff >= 0])]
				OUT[3] = median(all.durs[use.for.med])
				OUT[4:(3 + n.bands)] = all.mags
			}
		} else {
			OUT[1] = 4
		}
	}
	return(OUT)
}












