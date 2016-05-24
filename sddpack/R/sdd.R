sdd = function (A, kmax = 100, alphamin = 0.01, lmax = 100, rhomin = 10e-20) {
	betabar = 1
	yinit = 1 #only yinit=1 supported in this version
	idx = 1
	#initialize
	m = nrow(A)
	n = ncol(A)

	Xsav = matrix(0, nrow=m, ncol=0)
	Dsav = matrix(0, nrow=0, ncol=1)
	Ysav = matrix(0, nrow=n, ncol=0)
	rhosav = numeric(0)
	
	rho = sum(A^2)
	iitssav = numeric(0)
	#outer loop
	for (k in 1:kmax) {
		s = matrix(0, nrow=m, ncol=1)
		iits = 0
		while ((svd(s)$d[1])^2 < (rho / n)) {
			y = matrix(0, nrow=n, ncol=1)
			y[idx] = 1
			s = A %*% y
			if (k > 1) {
				s = s - (Xsav %*% (Dsav * (t(Ysav) %*% y)))
			}	
			idx = idx %% n + 1
			iits = iits + 1
		}
		iitssav[k] = iits	
				
		for (l in 1:lmax) {
			s = A %*% y
			if (k > 1) {
				s = s - (Xsav %*% (Dsav * (t(Ysav) %*% y)))
			}
	
			tmp = sddsolve(s, m)
			x = tmp[[1]]
			xcnt = tmp[[2]]
	
			s = t(A) %*% x
			if (k > 1) {
				s = s - (Ysav %*% (Dsav * (t(Xsav) %*% x)))
			}
	
			tmp2 = sddsolve(s, n)
			y = tmp2$x
			ycnt = tmp2$imax
			fmax = tmp2$fmax
	
			#check progress
			d = sqrt(fmax %*% ycnt) / (ycnt %*% xcnt)
		
			beta = d^2 %*% ycnt %*% xcnt
			if (l > 1) {
				alpha = (beta - betabar) / betabar
				if (alpha <= alphamin || abs(beta - betabar) < 10e-20) {
					break;
				}	
			}
			betabar = beta
		}
		
		#save
		Xsav = cbind(Xsav, x)
		Ysav = cbind(Ysav, y)
		Dsav = rbind(Dsav, d)
		rho = pmax(rho - beta, 0)
		rhosav[k] = rho

		if (rho <= rhomin) {
			break
		}
		if (k >= 2 &&  (rhosav[k] - rhosav[k-1] == 0)) {
			break
		}		
	}
	ret = list()
	colnames(Xsav)=rep("", ncol(Xsav))
	ret$x = Xsav
	ret$d = as.numeric(Dsav)
	colnames(Ysav)=rep("", ncol(Ysav))
	ret$y = Ysav
	return(ret)
}


sddsolve = function (s, m) {
	x = numeric(0)
	for (i in 1:m) {
		if (s[i] < 0) {
			x[i] = -1
			s[i] = -s[i]	
		} else {
			x[i] = 1	
		}
	}	

	sorts =-sort(-s)
	indexsort = order(-s)
	
	f = numeric(0)
	f[1] = sorts[1]
	for (i in 2:m) {
		f[i] = sorts[i] + f[i-1]	
	}
	f = f^2 / seq(1, m)
	imax = 1
	fmax = f[1]
	for (i in 2:m) {
		if (f[i] >= fmax) {
			imax = i
			fmax = f[i]	
		}
	}
	
	if ((imax + 1) <= m) {
		for (i in (imax+1):m) {
			x[indexsort[i]] = 0
		}
	}
	ret = list()
	ret$x = x
	ret$imax = imax
	ret$fmax = fmax
	return(ret)
}


