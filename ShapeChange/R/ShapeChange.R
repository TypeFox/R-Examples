changept = function(formula, family = gaussian(), data = NULL, k = NULL, knots = NULL, fir = FALSE, q = 3, pnt = FALSE, pen = 0, arp = FALSE, ci = FALSE, nloop = 1e+3, constr = TRUE, param = TRUE, gcv = FALSE)
{
	cl = match.call()
	if (is.character(family)) 
     		family = get(family, mode = "function", envir = parent.frame())
  	if (is.function(family)) 
     		family = family()
  	if (is.null(family$family)) 
     		stop("'family' not recognized!")
  	mf = match.call(expand.dots = FALSE)
  	m = match(c("formula", "data"), names(mf), 0L)
  	mf = mf[c(1L, m)]
  	mf[[1L]] = as.name("model.frame")
  	mf = eval(mf, parent.frame())
  	ynm = names(mf)[1]
  	mt = attr(mf, "terms")
    	y = model.response(mf, "any")
  	if (!is.null(names(y))) {
    		y = unname(y)
  	}
	if (family$family == "binomial") {
		if (class(y) == "factor") {
			y = ifelse(y == levels(y)[1], 0, 1)
		}
  	}
	if (family$family == "binomial" | family$family == "poisson") {
		wt.iter = TRUE
	} else {
		wt.iter = FALSE
	}
	zmat = NULL; zid = NULL; zid0 = NULL; zid1 = NULL; zid2 = NULL; znms = NULL; is_fac = NULL; vals = NULL; dist = 0
	for (i in 2:ncol(mf)) {
		if (is.character(attributes(mf[, i])$categ)) {
  			x = mf[, i]
  			categ = attributes(x)$categ
  			sh = attributes(x)$sh
			trd1 = attributes(x)$trd1
  			trd2 = attributes(x)$trd2
			up = attributes(x)$up
			xnm = attributes(x)$nm
		}
		if (is.null(attributes(mf[, i])$categ)) {
			if (!is.null(names(mf)[i])) {
				znms = c(znms, names(mf)[i])
				vals = c(vals, unique(mf[, i]))
			}
			#zmat = cbind(zmat, mf[, i])
			if (!is.matrix(mf[, i])) {
      				zid = c(zid, i)
        			if (is.factor(mf[, i])) {
	  				is_fac = c(is_fac, TRUE)
					ch_char = suppressWarnings(is.na(as.numeric(levels(mf[, i]))))
         			 	#if (any(ch_char)) {
          				#	vals = c(vals, unique(levels(mf[, i]))[2])
          				#} else {
	    				#	vals = c(vals, min(as.numeric(levels(mf[, i]))))
	  				#}
                 			nlvs = length(attributes(mf[, i])$levels)
         				zid0 = i + 0:(nlvs - 2) + dist
	  				zid1 = c(zid1, i + dist)
                 			zid2 = c(zid2, i + nlvs - 2 + dist)
       				   	dist = nlvs - 2
	  				#zmat0 = as.matrix(model.matrix(mt, mf)[, zid0], ncol = (length(zmat0) / length(y)))
					zmat0 = model.matrix(~ mf[, i])[, -1, drop = FALSE]
  				       	mat_cols = ncol(zmat0)
	  				mat_rm = NULL
  				       	rm_num = 0
	  				for (irm in 1:mat_cols) {
       	  					if (all(round(diff(zmat0[, irm]), 8) == 0)) {
                					mat_rm = c(mat_rm, irm)
          					}
   	  				}
	  				if (!is.null(mat_rm)) {
	  					zmat0 = zmat0[, -mat_rm, drop = FALSE]
                   				rm_num = rm_num + length(mat_rm)
				        }
	  				zmat = cbind(zmat, zmat0)
    			  	} else {
					is_fac = c(is_fac, FALSE)
					#zmat = cbind(zmat, mf[, i])
              				zmat = cbind(zmat, model.matrix(~ factor(mf[, i]))[, -1, drop = FALSE])			
					zid1 = c(zid1, i + dist)
					zid2 = c(zid2, i + dist)
            			}
			} else {
	 			is_fac = c(is_fac, FALSE)
  				zmat0 = mf[, i]
  				mat_cols = ncol(zmat0)
  				mat_rm = NULL
  				rm_num = 0
  				for (irm in 1:mat_cols) {
       	 				if (all(round(diff(zmat0[, irm]), 8) == 0)) {
               					mat_rm = c(mat_rm, irm)
         				}
  				}
  				if (!is.null(mat_rm)) {
 					zmat0 = zmat0[, -mat_rm, drop = FALSE]
					rm_num = rm_num + length(mat_rm)
  				}
  				zmat = cbind(zmat, zmat0)
  				#vals = c(vals, 1)
  				zid1 = c(zid1, i + dist)
  				zid2 = c(zid2, i + ncol(mf[, i]) - 1 + dist - rm_num)
  				zid = c(zid, i)
  				dist = ncol(mf[, i]) - 1
			}
		}
	}
#check
	if (!is.null(zmat)) {
		colnames(zmat) = paste(znms, vals[(length(vals) - ncol(zmat) + 1):length(vals)], sep = "")
	}
  	ans = changept.fit(x, y, zmat = zmat, family = family, categ = categ, m = k, knots = knots, sh = sh, fir = fir, q = q, pnt = pnt, pen = pen, arp = arp, ci = ci, nloop = nloop, trd1 = trd1, trd2 = trd2, up = up, constr = constr, wt.iter = wt.iter, param = param, gcv = gcv)
  	rslt = list(chpt = ans$chpt, knots = ans$knots, fhat = ans$fhat, fhat_x = ans$fhat_x, fhat_eta = ans$fhat_eta, fhat_eta_x = ans$fhat_eta_x, coefs = ans$bhat, zcoefs = ans$zcoefs, cibt = ans$cibt, categ = categ, sh = sh, x = x, y = y, xnm = xnm, znms = znms, ynm = ynm, m_i = ans$m_i, pos = ans$pos, sub = ans$sub, family = family, wt.iter = wt.iter, zmat = zmat, vals = vals, is_fac = is_fac, zid = zid, zid1 = zid1, zid2 = zid2, tms = mt, msbt = ans$msbt, bmat = ans$bmat, phi = ans$phi, sig = ans$sig, aics = ans$aics, lambda = ans$lambda, edf = ans$tr, edfu = ans$tru, lams = ans$lams, phisbt = ans$phisbt, sigsbt = ans$sigsbt)
  	rslt$call = cl
  	class(rslt) = "ShapeChange"
  	return (rslt) 
}

########
changept.fit = function(x, y, zmat = zmat, family = gaussian(), categ = categ, m = NULL, knots = NULL, sh = 1, fir = FALSE, q = 3, pnt = FALSE, pen = 0, arp = FALSE, lambdas = NULL, ci = FALSE, nloop = 1e+3, trd1 = -1, trd2 = -1, up = TRUE, constr = TRUE, wt.iter = wt.iter, param = TRUE, hs0 = NULL, ids0 = NULL, gcv = FALSE) {
	n = length(y)
	if (is.null(hs0) && arp) {  	
     		hs0 = 0:(n - 1)
     		id_1 = id_2 = hs0
     		id_mat = sapply(id_1, function(id_1) id_2 - id_1)
     		ids0 = sapply(hs0, function(h) {which(abs(id_mat) == h)})
	}
	if (categ == "inflect") {
		ans = ip.kts(x, y, zmat = zmat, m = m, knots = knots, q = q, pen = pen, pnt = pnt, arp = arp, lambdas = lambdas, pen_bt = NULL, p_bt = NULL, sh = sh, fir = fir, wt.iter = wt.iter, hs = hs0, ids = ids0, family = family, gcv = gcv) 
	} else if (categ == "mode") {
		ans = mode.kts(x, y, zmat = zmat, m = m, knots = knots, q = q, pen = pen, pnt = pnt, arp = arp, lambdas = lambdas, pen_bt = NULL, p_bt = NULL, sh = sh, wt.iter = wt.iter, hs = hs0, ids = ids0, family = family, gcv = gcv) 
	} else if (categ == "jp") {
		minsse = sum(y^2)
		for (i in 1:(n - 1)) {
			ansi = jp.pts(x, y, zmat = zmat, jpt = (sort(x))[i], i = NULL, m = m, knots = knots, q = q, pen = pen, pnt = pnt, up = up, trd1 = trd1, trd2 = trd2, constr = constr, wt.iter = wt.iter, arp = arp, lambdas = lambdas, pen_bt = NULL, p_bt = NULL, hs = hs0, ids = ids0, family = family, gcv = gcv) 		
			ssei = sum((y - ansi$fhat)^2)
			if (ssei < minsse) { 
				minsse = ssei; ans = ansi  
			}
		}
	} else {
		stop("change-point cannot be estimated!")
	}
	cibt = msbt = trs = trrs = lamsbt = sigsbt = psbt = NULL
	phisbt = list()
	if (ci & nloop > 0) {
		msbt = 1:nloop*0
		trs = trrs = lamsbt = sigsbt = psbt = 1:nloop*0
		phisbt = list()
		yvec = ans$fhat
#new
    		ord = order(x)
    		xs = sort(x)
    		x = xs
    		y = y[ord]
		evec = y - yvec
    		ys = matrix(0, nrow = nloop, ncol = n)
		ans_pen = ans_p = ans_lam = NULL
		for (j in 1:nloop) {
			if (family$family == "gaussian") {
           			if (!arp) {
				      esim = sample(evec, size = n, replace = TRUE)
           			} else {
             				ans_pen = ans$lambda
					ans_phi = ans$phi
					ans_sig = ans$sig
             				ans_p = length(ans$phi)
             				ans_lam = ans$lams
             				if (param) {
                				esim = as.vector(arima.sim(n = n, list(order = c(ans_p, 0, 0), ar = ans_phi), sd = sqrt(ans_sig)))
             				} else {
						#r0 = ans$r0
               					rmat = ans$rmat
						#rmat = rmat * r0
               					yvec = ans$fhat
               					evec = matrix(ans$es, ncol = 1)
               					umat = chol(rmat)
               					e_iid = as.vector(solve(t(umat), evec))
               					esamp = matrix(sample(e_iid, size = n, replace = TRUE), ncol = 1)
               					#esim = sqrt(ans_sig) * crossprod(umat, esamp)
             					esim = crossprod(umat, esamp)
						}
           			}
				ysim = yvec + esim
			} else if (family$family == "binomial") {
				ysim  = rbinom(n, size = 1, prob = yvec)
			} else if (family$family == "poisson") {
				ysim  = rpois(n, lambda = yvec)
			} 
			if (categ == "inflect") {
				ansj = try(ip.kts(x, ysim, zmat = zmat, m = m, knots = knots, q = q, pen = pen, pnt = pnt, arp = arp, lambdas = ans_lam, pen_bt = ans_pen, p_bt = ans_p, sh = sh, fir = fir, wt.iter = wt.iter, hs = hs0, ids = ids0, family = family, gcv = gcv))
        			if (class(ansj) == "try-error") next 
			} else if (categ == "mode") {            
				ansj = try(mode.kts(x, ysim, zmat = zmat, m = m, knots = knots, q = q, pen = pen, pnt = pnt, arp = arp, lambdas = ans_lam, pen_bt = ans_pen, p_bt = ans_p, sh = sh, wt.iter = wt.iter, hs = hs0, ids = ids0, family = family, gcv = gcv))
        			if (class(ansj) == "try-error") next 
			} else if (categ == "jp") {
				minsse = sum(ysim^2)				
				for (i in 1:(n - 1)) {
					ansi = try(jp.pts(x, ysim, zmat = zmat, jpt = sort(x)[i], i = NULL, m = m, knots = knots, q = q, pen = pen, pnt = pnt, up = up, trd1 = trd1, trd2 = trd2, constr = constr, wt.iter = wt.iter, arp = arp, lambdas = lambdas, pen_bt = ans_pen, p_bt = ans_p, hs = hs0, ids = ids0, family = family, gcv = gcv))
					if (class(ansi) == "try-error") next 
#new:
					if (arp) {
						rmat = ansi$rmat 
						ssei = crossprod((ysim - ansi$fhat), chol2inv(chol(rmat))) %*% (ysim - ansi$fhat)
					} else {
						ssei = sum((ysim - ansi$fhat)^2)
					}
					if (ssei < minsse) {
						minsse = ssei; ansj = ansi  
					}
				}
			}
			msbt[j] = ansj$chpt
			if (arp) {
				sigsbt[j] = ansj$sig
				phisbt[[j]] = ansj$phi
				lp = length(ansj$phi)
				if (lp > 1) {
					psbt[j] = max(which(round(ansj$phi, 6) != 0)) 
				} else {
					if (round(ansj$phi, 6) != 0) {
						psbt[j] = 1
					} else {psbt[j] = 0}
				}    						
				trs[j] = ansj$tr
				trrs[j] = ansj$trr
      				lamsbt[j] = ansj$lambda
			}
  	}
	cibt = quantile(sort(msbt), probs = c(.025, .975))
	}
	rslt = list(chpt = ans$chpt, knots = ans$knots, fhat = ans$fhat, fhat_x = ans$fhat_x, fhat_eta = ans$fhat_eta, fhat_eta_x = ans$fhat_eta_x, sse = ans$sse, df = ans$df, pos = ans$pos, dist = ans$dist, bhat = ans$bhat, zcoefs = ans$zcoefs, pvals.beta = ans$pvals.beta, se.beta = ans$se.beta, tval = ans$tval, cibt = cibt, m_i = ans$m_i, sub = ans$sub, tr = ans$tr, trr = ans$trr, tru = ans$tru, msbt = msbt, phisbt = phisbt, psbt = psbt, sigsbt = sigsbt, bmat = ans$bmat, dmat = ans$dmat, phi = ans$phi, sig = ans$sig, aics = ans$aics, aiclst = ans$aiclst, es = ans$es, lambda = ans$lambda, lams = ans$lams, edfs = ans$edfs, gcvus = ans$gcvus, lambdas_pen = ans$lambdas_pen, rmat = ans$rmat, aicmat = ans$aicmat, sigmat = ans$sigmat, trsmat = ans$trsmat, trrsmat = ans$trrsmat, psmat = ans$psmat, fsmat = ans$fsmat, phismat = ans$phismat, psmat = ans$psmat, pensmat = ans$pensmat, hs = hs0, ids = ids0)
	rslt
}

######################
jp = function(x, trd1 = -1, trd2 = -1, up = TRUE) {
  cl = match.call()
  pars = match.call()[-1]
  attr(x, "nm") = deparse(pars$x)
  attr(x, "categ") = "jp"
  attr(x, "trd1") = trd1
  attr(x, "trd2") = trd2
  attr(x, "up") = up
  x
}

tp = function(x, sh = 1) {
  cl = match.call()
  pars = match.call()[-1]
  attr(x, "nm") = deparse(pars$x)
  attr(x, "categ") = "mode"
  attr(x, "sh") = sh
  x
}

ip = function(x, sh = 1) {
  cl = match.call()
  pars = match.call()[-1]
  attr(x, "nm") = deparse(pars$x)
  attr(x, "categ") = "inflect"
  attr(x, "sh") = sh
  x
}

######################
fitted.ShapeChange = function(object, ...) {
  if (!inherits(object, "ShapeChange")) { 
    warning("calling fitted.ShapeChange(<fake-ShapeChange-object>) ...")
  }
  ans = as.vector(object$fhat)
  ans
}

coef.ShapeChange = function(object, ...) {
  if (!inherits(object, "ShapeChange")) { 
    warning("calling coef.ShapeChange(<fake-ShapeChange-object>) ...")
  }
  ans = object$coefs
  ans	
}

confint.ShapeChange = function(object, parm = NULL, level = NULL, ...) {
  if (!inherits(object, "ShapeChange")) { 
    warning("calling confint.ShapeChange(<fake-ShapeChange-object>) ...")
  }
  ans = object$cibt
  ans	
}

residuals.ShapeChange = function(object, ...) {
  if (!inherits(object, "ShapeChange")) { 
    warning("calling resid.ShapeChange(<fake-ShapeChange-object>) ...")
  }
  ans = object$y - object$fhat
  ans	
}

##########################
plot.ShapeChange = function(x, y = NULL, tt = TRUE, xlab = NULL, ylab = NULL, color = "mediumorchid4", type = "p", lty = 1, lwd = 2, cex = .5, has.leg = FALSE, loc = NULL, rugged = FALSE, detailed = FALSE,...) {
	if (!inherits(x, "ShapeChange")) { 
		warning("calling plot.ShapeChange(<fake-ShapeChange-object>) ...")
 	}
 	object = x
 	x = object$x
	zmat = object$zmat
	znm = object$znm
	zid = object$zid
	zid1 = object$zid1 
	zid2 = object$zid2
	zbhat = object$zcoefs
  	tms = object$tms
	vals = object$vals
	znm_2 = dimnames(zmat)[[2]]
	xnm = object$xnm
	ord = order(x)
	x0s = sort(x)
	x = (x0s - min(x0s)) / (max(x0s) - min(x0s))
	y = object$y
	ynm = object$ynm
	y = y[ord]
	sub = object$sub
	fhat = object$fhat
	fhat_x = object$fhat_x
	fhat_eta_x = object$fhat_eta_x
	is_fac = object$is_fac
	bhat = object$coefs
  	chpt = object$chpt	
	pos = object$pos
#constr = object$constr
	ci = object$ci
	m_i = object$m_i
  	categ = object$categ
	sh = object$sh
	knots = object$knots
	knots = (max(x0s) - min(x0s)) * knots + min(x0s)
	k = length(knots)
	lb = length(bhat)
	wt.iter = object$wt.iter
	family = object$family
#print (family$family)
	cicfamily = CicFamily(family)
	muhat.fun = cicfamily$muhat.fun
	if (!is.null(zmat)) {
		zu = sort(unique(zmat))
	}
	palette = c("lightblue", "peachpuff", "limegreen", "grey", "wheat", "yellowgreen", "seagreen1", "palegreen", "azure", "whitesmoke")
  	ylim = c(min(c(min(y), min(fhat))), max(c(max(y), max(fhat))))
#new:
	if (categ == "mode" | categ == "jp") {
		ans0 = bqspl(x, k, pic = TRUE, spl = FALSE)
	} else if (categ == "inflect") {
		ans0 = bcspl(x, k, pic = TRUE, spl = FALSE)
	}
	xpl = ans0$xpl
	xpl0 = 0:1000/1000*(max(x0s) - min(x0s)) + min(x0s)
	bpl = ans0$bpl    
	if (categ == "jp") {
		n = length(x)
		ijp = pos
		d1 = d2 = 1:1001*0
		d1[xpl <= (x[ijp] + x[ijp + 1]) / 2] = (2 * ijp + 1) / 2 - n
		d1[xpl > (x[ijp] + x[ijp + 1]) / 2] = (2 * ijp + 1) / 2
		d2[xpl <= (x[ijp] + x[ijp + 1]) / 2] = 0
		d2[xpl > (x[ijp] + x[ijp + 1]) / 2] = xpl[xpl > (x[ijp] + x[ijp+1]) / 2] - (x[ijp] + x[ijp+1]) / 2
		d2 = d2 - m_i
		if (ijp == 1 | ijp == (n - 1)) {
			d2 = NULL
		}
		di = cbind(d1, d2)	
		bpl = cbind(bpl, di)
	} 
	if (is.null(zmat)) {
		thp = bpl %*% bhat
	} else {
		thp = bpl %*% as.vector(bhat)[1:(ncol(bpl))]	
	}	
	if (tt) {
    		if (categ == "mode") {
     	 		tt0 = "Mode Estimate: "
    	} else if (categ == "jp") {
      		tt0 = "Jump-Point Estimate: "
    	} else if  (categ == "inflect") {
      		tt0 = "Inflection-Point Estimate: "
    	} else {
      		stop ("Wrong Shape!!")
    	}
    		main = paste(tt0, round(chpt, 2))
  	} else { 
    		main = NULL
  	}
	thp0 = thp
	if (wt.iter & categ != "inflect") {
		thp = muhat.fun(thp0, fml = family$family)
	}
	r1 = range(y); r2 = range(thp); r3 = range(fhat)
	rmin = min(c(r1, r2, r3)); rmax = max(c(r1, r2, r3))
	rg = rmax - rmin 
  	if (!is.null(xlab)) {
     		xnm = xlab
  	}
  	if (!is.null(ylab)) {
     		ynm = ylab
  	}
	plot(x0s, y, type = "n", col = 1, cex = cex, ylab = ynm, xlab = xnm, ylim = c(rmin, rmax), main = main)
	if (family$family == "binomial" & rugged) {
		#loc = "right"
		rug(x0s[y == 0])
		rug(x0s[y == 1], side = 3)
	} else {
		points(x0s, y, lwd = 1, cex = cex)
	}	
	if (is.null(zmat)) {
		if (detailed) {
			lines(xpl0, thp, type = 'l', lty = 2, xlab = '', ylab = '', col = color, lwd = 2)
			points(x0s, fhat, col = color, cex = cex)   	
		} else {
  			lines(x0s, fhat, type = 'l', lty = 2, xlab = '', ylab = '', col = color, lwd = lwd)
		}		
	} else {
		lzb = length(zbhat)
		if (detailed) {
			lines(xpl0, thp, type = 'l', lty = 2, xlab = '', ylab = '', col = color, lwd = 2)	
			points(x0s, fhat_x, col = color, cex = cex)
			for (i in 1:lzb) {                                      
				thpi = thp + zbhat[i]
				if (wt.iter & categ != "inflect") {
					thpi = muhat.fun(thpi, fml = family$family)
					fhat_eta_xi = fhat_eta_x + zbhat[i]
					fhat_xi = muhat.fun(fhat_eta_xi, fml = family$family)				
				} else {
					fhat_xi = fhat_x + zbhat[i]
				}
				lines(xpl0, thpi, type = 'l', lty = 2, xlab = '', ylab = '', col = palette[i], lwd = 2)	
				points(x0s, fhat_xi, col = palette[i], cex = cex)			
			}
		} else {
  			lines(x0s, fhat_x, type = 'l', lty = 2, col = color, lwd = 2)
			for (i in 1:lzb) {                                      
				if (wt.iter & categ != "inflect") {
					fhat_eta_xi = fhat_eta_x + zbhat[i]
					fhat_xi = muhat.fun(fhat_eta_xi, fml = family$family)				
				} else {
					fhat_xi = fhat_x + zbhat[i]
				}
				lines(x0s, fhat_xi, type = 'l', lty = 2, col = palette[i], lwd = 2)				
			}
		}
	}
	#text(chpt, ylim[1], round(chpt, 4), col = color, lwd = 2)
	if (!is.null(ci)) {
		lvl = 95
		if (is.null(loc)) {
			if (wt.iter) {
				loc = "left"
			} else {
				loc = "topright"
			} 
		}
		legend(loc, bty = "n", paste(lvl, "%", c("Bootstrap C.I.  "), "[", round(ci[1], 4), ",", round(ci[2], 4), "]"), col = color, cex = 1)
		abline(v = ci[1], lty = 2)
		abline(v = ci[2], lty = 2)
	}
	if (has.leg) {
		if (is.null(zmat)) {
  			legend("topleft", bty = "n", c("f_hat"), col = color, lty = 2, lwd = 2)
		} else {
			legend(x = min(x0s), y = rmax , bty = "n", c("constrained f_hat"), col = color, lty = 2, lwd = 2)
			#if (!is_fac) {
			#	znm_3 = znm_2
			#} else {	
			#	znm_3 = paste("covariate the ", 1:lzb, "th column")
			#}
			#if (lzb >= 1) {
			#	for (i in 1:lzb) {	
			#		legend(x = min(x0s), y = rmax - i * .12 * rmax, bty = "n", znm_3[i], col = palette[i], lty = 1, lwd = 2)
			#	}
			#}
		}
  	}
  	#points(knots, 1:length(knots)*0 + rmin, pch = 'X')
  	#if (categ == "jp") {
        #	abline(v = sub[1], lty = 2, col = 2, lwd = 2)
    	#	abline(v = sub[2], lty = 2, col = 2, lwd = 2)
  	#}
}


###########################
#print.summary.ShapeChange#
###########################
print.summary.ShapeChange = function(x,...) {
	if (!is.null(x$zcoefs)) {
		cat("Call:\n")
		print(x$call)
		cat("\n")
		cat("Coefficients:")
		cat("\n")
		printCoefmat(x$coefficients, P.values = TRUE, has.Pvalue = TRUE)
		cat("\n")
	} else {
		print ("No predictor is defined")	
	}
}

#####################
#summary.ShapeChange#
#####################
#summary.ShapeChange = function(object,...) {
#	if (!is.null(object$zcoefs)) {
#		family = object$family 
#		wt.iter = object$wt.iter
#		coefs = object$zcoefs
#		se = object$se.beta
#		tval = object$tval
#    		pvalbeta = object$pvals.beta
#		n = length(coefs)
#		zmat = object$zmat
#		is_fac = object$is_fac
#		if (wt.iter) {
#			rslt1 = data.frame("Estimate" = round(coefs, 4), "StdErr" = round(se, 4), "z.value" = round(tval, 4), "p.value" = round(pvalbeta, 4))			
			#rownames(rslt1)[1] = "(Intercept)"
#		} else {
#			rslt1 = data.frame("Estimate" = round(coefs, 4), "StdErr" = round(se, 4), "t.value" = round(tval, 4), "p.value" = round(pvalbeta, 4))
#		}
#		rownames(rslt1)[1:n] = colnames(zmat)
#		rslt1 = as.matrix(rslt1)
#		ans = list(call = object$call, coefficients = rslt1, zcoefs = coefs, family = family)
		#class(ans) = "summary.ShapeChange"
		#ans
#	} else {
#		ans = list(zcoefs = object$zcoefs)
		#class(ans) = "summary.ShapeChange"
		#ans
#	}
#	class(ans) = "summary.ShapeChange"
#	ans
#}

#####################
#predict.ShapeChange#
#####################
predict.ShapeChange = function(object, newData,...) {
	family = object$family
	cicfamily = CicFamily(family)
	muhat.fun = cicfamily$muhat.fun	
	coefs = object$coefs
	zcoefs = object$zcoefs
  	categ = object$categ
	sh = object$sh
  	tt = object$tms
 	knots = object$knots
  	nk = length(knots)
  	x = object$x
  	xs = sort(x)
  	xs = (xs - min(xs)) / (max(xs) - min(xs))
  	wt.iter = object$wt.iter
  	bmat = object$bmat
  	zmat = object$zmat
  	vals = object$vals
  	add = FALSE
	if (!inherits(object, "ShapeChange")) { 
    		warning("calling predict.ShapeChange(<fake-ShapeChange-object>) ...")
  	}
	if (missing(newData) | is.null(newData)) {
		ans = list(fhat = object$fhat, fhat_x = object$fhat_x)
		return (ans) 
	}
	Terms = delete.response(tt)
	m = model.frame(Terms, newData)
	newdata = m
  	nx = newdata[, 1]
  	if (!is.null(zcoefs)) {
     		nz = newdata[, 2]
     		if (!any(nz %in% vals)) {
        		stop ('New categorical variable value not found in the fit!')
     		}
#check!
     		if (length(unique(nz)) > 1) {
		        zmat = model.matrix(~ factor(nz))[, -1, drop = FALSE]
	     	} else {
	    		pos = which(vals == unique(nz))
       			add = TRUE
     		}
#print (head(zmat))
	}
  	nx0 = nx
  	snx = sort(nx)
#print (nx)
	nx = (snx - min(x)) / (max(x) - min(x))
  	if (any(nx0 > max(x)) | any(nx0 < min(x))) {
     		warning ('User is extrapolating!')
  	} else {
      		nx = (snx - min(x)) / (max(x) - min(x))
  	}
  	if (categ == 'mode' | categ == 'jp') {
     		ans0 = bqspl(x = nx, m = nk, knots = knots, pic = FALSE, x0 = xs)
  	} else if (categ == 'inflect') {
    		ans0 = bcspl(x = nx, m = nk, knots = knots, pic = FALSE, x0 = xs)
  	}
  	nbmat = ans0$bmat
	if (is.null(zcoefs)) {
     		nthb = nbmat %*% coefs
     		nthb_x = nthb
  	} else {
    		if (categ == 'mode' | categ == 'jp') {
       			nthb_x = nbmat[, 1:(nk + 1), drop = FALSE] %*% coefs[1:(nk + 1)]
    		} else if (categ == 'inflect') {
       			nthb_x = nbmat[, 1:(nk + 2), drop = FALSE] %*% coefs[1:(nk + 2)]
    		}
    		if (!add) {
       			nbmat = cbind(nbmat, zmat)
       			nthb = nbmat %*% coefs
    		} else {
      			if (pos == 1) {nthb = nthb_x} else {nthb = nthb_x + zcoefs[pos - 1]}
    		}
	}   
	if (wt.iter & categ == 'mode') {
     		if (is.null(zcoefs)) {
        		nthb = muhat.fun(nthb, fml = family$family)
        		nthb_x = nthb
     		} else {
       			nthb_x = muhat.fun(nthb_x, fml = family$family)
       			nthb = muhat.fun(nthb, fml = family$family)
     		}
  	}
	ans = list(fhat = nthb, fhat_x = nthb_x, nbmat = nbmat)
	return (ans) 
 }

############################
#find penalty terms #
############################
find_pen = function(aims = NULL, Q = NULL, B = NULL, D = NULL, PNT = TRUE, Y = NULL, D0 = NULL, GCV = FALSE) {
        la = length(aims)
        lam_use = lambdas = NULL
	gcvus = NULL
	if (!GCV) {         
		if (PNT) {
        	    for (i in 1:la) {
        	      x.l = 1e+10; x.r = 1e-10                        
        	      f = function(pen0) sum(diag(B %*% solve((Q + pen0 * D), t(B)))) - aims[i]
        	      if (abs(f(x.r)) < 1e-4) {
        	         lambda = 0
        	      } else {
        	        lambda = uniroot(f, c(x.l, x.r), tol = 1e-8)$root
        	      }
        	      lambdas = c(lambdas, lambda)
        	   }
        	 } else {
        	   lambdas = 0
		 }
         	 #lambdas
	} else {
		n = length(Y)
		for (i in 1:la) {
#x.r is the smallest possible lambda
        	 	x.l = 1e+10; x.r = 1e-10                        
        		f = function(pen0) sum(diag(B %*% solve((Q + pen0 * D), t(B)))) - aims[i]
        		if (abs(f(x.r)) < 1e-4) {
        	        	lambda = 0
        	      	} else {
        	        	lambda = uniroot(f, c(x.l, x.r), tol = 1e-8)$root
        	      	}
			pmatu = B %*% solve((crossprod(B) + lambda * crossprod(D0)), t(B)) 
		 	muhatu = pmatu %*% Y
		 	sseu = sum((Y - muhatu)^2)
		 	#edfu = sum(diag(pmat))
			edfu = aims[i]
		 	gcvu = sseu / (1 - edfu/n)^2
        	      	lambdas = c(lambdas, lambda)
			gcvus = c(gcvus, gcvu)
        	 }
	}
	if (!GCV) {
		lam_use = lambdas
	} else {
		lam_use = lambdas[gcvus == min(gcvus)]
		if (length(lam_use) > 1) {
			stop ('Check find_pen!')
		}
	}
        rslt = list(lam_use = lam_use, lambdas = lambdas, gcvus = gcvus)
	rslt 
}


#############
#mode.kts#
#############
mode.kts = function(x, y, zmat = NULL, m = NULL, knots = NULL, sh = 1, q = 3, pen = 0, pnt = FALSE, wt.iter = FALSE, arp = FALSE, lambdas = NULL, pen_bt = NULL, p_bt = NULL, hs = NULL, ids = NULL, family = gaussian(), gcv = FALSE) {
	linkfun = family$linkfun
	cicfamily = CicFamily(family)
	llh.fun = cicfamily$llh.fun
	etahat.fun = cicfamily$etahat.fun
	gr.fun = cicfamily$gr.fun
	wt.fun = cicfamily$wt.fun
	zvec.fun = cicfamily$zvec.fun
	muhat.fun = cicfamily$muhat.fun
	ysim.fun = cicfamily$ysim.fun
	deriv.fun = cicfamily$deriv.fun
	dev.fun = cicfamily$dev.fun 
	xs = sort(x)	
	ord = order(x)
  	#if (!no.ord) {
  	y = y[ord]
  	#}
	if (sh == -1 & !wt.iter) {
		y = -y
	}	
	if (is.matrix(y)) {
		y = as.numeric(y)
	}
#new:
	x = (xs - min(xs)) / (max(xs) - min(xs))
	n = length(x)	
	sm = 1e-5
	if (is.null(m)) {
		m0 = 3 + round(n^(1 / 7)) 
		m1 = m0
		m = NULL
		while (is.null(m)) { 
			pts = floor((n - 2) / (m1 - 1))
			rem_pts = (n - 2) %% (m1 - 1)
			if (pts > 2) {
				m = m1
			} else if (pts == 2) {
				if (rem_pts / (m1 - 1) >= 1) {
					m = m1
				} else {
					m1 = m1 - 1	
				}
			} else {
				m1 = m1 - 1
			}
		}
	}
	m2 = m
	ans = bqspl(x, m = m2, knots = knots, pnt = pnt)
	knots = ans$knots
	bmat = ans$bmat
	#nkts = length(knots)
#new
	if (!is.null(zmat)) {
		bmat = cbind(bmat, zmat)
    		np = ncol(zmat)
	} else {np = 0}
  	qv0 = crossprod(bmat)
  	qv = qv0
  	slopes = ans$slopes
	m = length(bmat) / n
  	la = 1
#new: return edfs
	edfs = NULL
	gcvus = NULL
	lambda = lambdas = lambdas_pen = NULL
	if (pnt) {
		dmat = matrix(0, nrow = (m - q), ncol = m)
#  third-order
		if (q == 3) {
			for (i in 4:m) {
				dmat[i-3, i-3] = 1; dmat[i-3, i-2] = -3; dmat[i-3, i-1] = 3; dmat[i-3, i] = -1
			}
		}
# second order
		if (q == 2) {
			for (i in 3:m) {
				dmat[i-2, i-2] = 1; dmat[i-2, i-1] = -2; dmat[i-2, i] = 1
			}
		}
# first order
		if (q == 1) {
			for (i in 2:m) {
				dmat[i-1, i-1] = 1; dmat[i-1, i] = -1
			}
		}
# zero order
		if (q == 0) {
			for (i in 1:m) {
				dmat[i, i] = 1
			}
		}
    		if (!wt.iter) {
       			dv0 = crossprod(dmat)
       			if (!arp) {
          			#if (pen == 0) {
             			#	pen = find_pen(aims = 6, Q = qv0, B = bmat, D = dv0, PNT = pnt)
          			#}
				if (pen == 0) {
					if (m < 6) {
						if (!gcv) {
# the unconstr number of columns of bmat 
							edfs = m
						} else {
							edfs = m:(m+2)
						}
						ans_pen = find_pen(aims = edfs, Q = qv0, B = bmat, D = dv0, PNT = pnt, Y = y, D0 = dmat, GCV = gcv)
	             				pen = ans_pen$lam_use
						gcvus = ans_pen$gcvus
						lambdas_pen = ans_pen$lambdas 
						
					} else {
						if (!gcv) {
# the unconstr number of columns of bmat 
							edfs = 6
						} else {
							edfs = 6:8
						}
						ans_pen = find_pen(aims = edfs, Q = qv0, B = bmat, D = dv0, PNT = pnt, Y = y, D0 = dmat, GCV = gcv)
						pen = ans_pen$lam_use
						gcvus = ans_pen$gcvus
						lambdas_pen = ans_pen$lambdas 
					}
		        	}
          			lambdas = pen 
       			}
       			if (arp & is.null(p_bt)) {        
	        		if (m < 10) {
            				aims = 5:m
         			} else {
            				aims = 5:10
         			}
				edfs = aims 
         			if (is.null(lambdas)) {
            				if (pen == 0) {
						ans_pen = find_pen(aims = aims, Q = qv0, B = bmat, D = dv0, PNT = pnt, Y = y, D0 = dmat, GCV = gcv)
               					lambdas = ans_pen$lam_use
            				} else {
               					lambdas = pen
            				}
         			}
       			} 
      			if (!is.null(p_bt)) {
          			lambdas = as.vector(pen_bt)
       			}
       			la = length(lambdas)
    		}
	} else {
		dmat = NULL 
	}
	dfmat = bmat %*% solve(qv, t(bmat))
  	pos = pos2 = NULL
	tr = trr = sig = aics = es = rmat = phi = minsse = NULL
	aiclst = aicmat = sigmat = trsmat = trrsmat = trslst = psmat = fsmat = phismat = psmat = pensmat = NULL
  	tr_use0 = tr_use = trr_use = NULL
	if (!wt.iter) {
		cv = crossprod(bmat, y)
    		smat = slopes
		if (!is.null(zmat)) {
			nilmat = matrix(0, nrow = nrow(smat), ncol = ncol(zmat))
			smat = cbind(smat, nilmat)
		}
#search for a window 
    		for (ila in 1:la) {
        		if (pnt) {
           			pen = lambdas[ila]
           			qv = qv0 + pen * dv0
        			dfmat = bmat %*% solve(qv, t(bmat))
				tr_use0 = c(tr_use0, sum(diag(dfmat)))
        		} else {qv = qv0}   
        		minsse = sum(y^2)
        		for (j in 1:m2) {
          			sl = smat 
          			sl[j:m2, ] = -sl[j:m2, ]         
          			qans = qprog(qv, cv, sl, 1:nrow(sl)*0, msg = FALSE)
          			bh = fitted(qans)
          			fhat = bmat %*% bh
          			fhat_x = bmat[, 1:(m2 + 1), drop = FALSE] %*% bh[1:(m2 + 1), , drop = FALSE]
				sse = mean((y - fhat)^2)
				if (sse < minsse) {thb = fhat; thb_x = fhat_x; minsse = sse; pos = j; bhat = bh; tr = tr_use0; lambda = pen}
       			} 
       			pos2 = c(pos2, pos)
       			pos_use = round(mean(pos2))
		} 
    		if (arp) {
        		st = pos_use - 2
        		ed = pos_use + 2   
			#st = 1
			#ed = m2        		
			if (st < 1) st = 1
        		if (ed > m2) ed = m2
			#minllh = 1e+100 #sum(y^2)
			minaic = sum(y^2)			
			aiclst = fslst = siglst = phislst = penslst = trslst = trrslst = list()
			iterlst = 0
			miniter = NULL
        		for (j in st:ed) {
				iterlst = iterlst + 1
            			sl = smat 
		      		sl[j:m2, ] = -sl[j:m2, ]
            			qans = make_arp(y, pnt = pnt, dmat = dmat, bmat = bmat, sl = sl, p_max = 2, m2 = m2, lambdas = lambdas, pen_fit = pen_bt, p_fit = p_bt, dv0 = dv0, qv0 = qv0, cv = cv, hs = hs, ids = ids)
				bh = qans$thetahat
            			theta = qans$fhat_use
				theta_x = bmat[, 1:(m2 + 1), drop = FALSE] %*% bh[1:(m2 + 1), , drop = FALSE] 
				aiclst[[iterlst]] = qans$aics	
				siglst[[iterlst]] = qans$sigs
				trslst[[iterlst]] = qans$trs
				trrslst[[iterlst]] = qans$trrs
				fslst[[iterlst]] = qans$fhats
				phislst[[iterlst]] = qans$phis
				penslst[[iterlst]] = qans$pens 
            			r_use = qans$r_use
				phi_use = qans$phi_use
            			sig_use = qans$sig_use
            			aics_use = qans$aics
            			tr_use = qans$tr_use
				trr_use = qans$trr_use
            			es_use = qans$es_use
            			lambda_use = qans$lambda_use       								
				aicj = qans$aicmin
	          		if (aicj < minaic) {thb = theta; thb_x = theta_x; minaic = aicj; pos = j; bhat = bh; miniter = iterlst; phi = phi_use; sig = sig_use; aics = aics_use; tr = tr_use; trr = trr_use; es = es_use; lambda = lambda_use; rmat = r_use}
        		}         
			aicmat = aiclst[[miniter]]
			sigmat = siglst[[miniter]] 
			trsmat = trslst[[miniter]]
			trrsmat = trrslst[[miniter]]  
			fsmat = fslst[[miniter]]   
			phismat = phislst[[miniter]] 
			pensmat = penslst[[miniter]]  
			id_minaic = which(aicmat == min(aicmat))
			aics = aicmat 
			#sig_tst = sigmat[id_minaic]
			#tr_tst = trsmat[id_minaic]
			#lambda_tst = pensmat[id_minaic] 			                         		
		}
	} else {
		if (!is.null(zmat)) {
			nc = m2 + 1 + ncol(zmat)
       		} else {nc = m2 + 1}
		smat = slopes
		if (!is.null(zmat)) {
			nilmat = matrix(0, nrow = nrow(smat), ncol = ncol(zmat))
			smat = cbind(smat, nilmat)
       		}
       		qv0 = crossprod(bmat)
       		qv = qv0
       		if (pnt) {       
#print (pnt)    
          		dv0 = crossprod(dmat)
          		if (pen == 0) {
				if (m < 6) {
					if (!gcv) {
# the unconstr number of columns of bmat 
						edfs = m
					} else {
						edfs = m:(m+2)
					}
					ans_pen = find_pen(aims = edfs, Q = qv0, B = bmat, D = dv0, PNT = pnt, Y = y, D0 = dmat, GCV = gcv)
             				pen = ans_pen$lam_use
					gcvus = ans_pen$gcvus
					lambdas_pen = ans_pen$lambdas 
				} else {
					if (!gcv) {
# the unconstr number of columns of bmat 
						edfs = 6
					} else {
						edfs = 6:8
					}
					ans_pen = find_pen(aims = edfs, Q = qv0, B = bmat, D = dv0, PNT = pnt, Y = y, D0 = dmat, GCV = gcv)
					pen = ans_pen$lam_use
					gcvus = ans_pen$gcvus
					lambdas_pen = ans_pen$lambdas 
				}
	        	}
			lambda = pen; lambdas = pen 
          		qv = qv0 + pen * crossprod(dmat)
#print (pen)
       		}
       		dfmat = bmat %*% solve(qv, t(bmat))
       		tr = sum(diag(dfmat))
		minllh = 1e+3
		for (j in 1:m2) {
			etahat = etahat.fun(n, y, fml = family$family)
			gr = gr.fun(y, etahat, weights = NULL, fml = family$family)  
			wt = wt.fun(etahat, n, weights = NULL, fml = family$family)
			cvec = wt * etahat - gr
			qvk = crossprod(bmat, diag(wt)) %*% bmat
           		if (pnt) {
				qvk = qvk + pen * dv0
           		}
			cveck = crossprod(bmat, cvec)			
			sl = smat 
			sl[j:m2, ] = -sl[j:m2, ]
			qans = qprog(qvk, cveck, sl, 1:m2*0, msg = FALSE)
			bh = qans$thetahat
			etahat = bmat %*% bh
			muhat = muhat.fun(etahat, fml = family$family)
			diff = 1
			nrep = 0
			while (diff > sm & nrep < 100) {
				oldmu = muhat	
				nrep = nrep + 1
				gr = gr.fun(y, etahat, weights = NULL, fml = family$family)	
				wt = wt.fun(etahat, n, weights = NULL, fml = family$family) 
				cvec = wt * etahat - gr
				qvk = crossprod(bmat, diag(wt)) %*% bmat
				if (pnt) {
					qvk = qvk + pen * dv0
#print (pen)
		            	}
				cveck = crossprod(bmat, cvec)
				qans = qprog(qvk, cveck, sl, 1:m2*0, msg = FALSE)
				bh = qans$thetahat
				etahat = bmat %*% bh
				muhat = muhat.fun(etahat, fml = family$family)
				diff = mean((muhat - oldmu)^2)	
         		}
			etahat_x = bmat[, 1:(m2 + 1), drop = FALSE] %*% bh[1:(m2 + 1), ,drop = FALSE]
			muhat_x = muhat.fun(etahat_x, fml = family$family)
			d_use = t(qans$xmat)
			qdf = qans$df
			tdf = ncol(bmat)
			imat = diag(tdf)
			if (qdf < tdf) {
   				pd = d_use %*% solve(crossprod(d_use), t(d_use))
   				pmat0 = imat - pd
			} else {
  				pmat0 = imat
			}   
			umat = chol(qvk)
			iumat = diag(ncol(umat))
			uinv = backsolve(umat, iumat)
			bu = bmat %*% uinv
			pmat = bu %*% tcrossprod(pmat0, bu) %*% diag(wt)
			tr = sum(diag(pmat))
         		llh = llh.fun(y, muhat, etahat, n, weights = NULL, fml = family$family)
			if (llh < minllh) {thb = muhat; thb_x = muhat_x; thb_eta = etahat; thb_eta_x = etahat_x; minllh = llh; pos = j; bhat = bh; tr = tr}
		}
	}
  	dist = round(slopes %*% bhat[1:(m2 + 1)], 6)
	sub = NULL 
	if (abs(dist[pos]) < sm) {
		if (pos == 1) {
			obs.r = 2:m2
			dist.r = dist[obs.r]
			if (dist.r[1] < (-sm)) {
				x.l = knots[pos]; x.r = knots[pos]
			} else if (abs(dist.r[1]) < sm) {
				if (any(dist.r < (-sm))) {
					id.r = (obs.r[dist.r < (-sm)])[1] - 1
					x.l = knots[pos]; x.r = knots[id.r]
				} else {
					x.l = knots[pos]; x.r = knots[m2]
				}
			}
		} else if (pos == m2) {
			obs.l = 1:(m2 - 1)
			dist.l = dist[obs.l]
			if (dist.l[m2 - 1] > sm) {
				x.l = knots[pos]; x.r = knots[pos]
			} else if (abs(dist.l[m2 - 1]) < sm) {
				if (any(dist.l > sm)) {
					id.l = tail(obs.l[dist.l > sm], 1) + 1
					x.l = knots[id.l]; x.r = knots[pos]
				} else {
					x.l = knots[1]; x.r = knots[pos]
				}
			}
		} else {
			obs.l = 1:(pos - 1); obs.r = (pos + 1):m2
			dist.l = dist[obs.l]; dist.r = dist[obs.r]
			if (dist.l[pos - 1] > sm) {
				x.l = knots[pos]
			} else if (abs(dist.l[pos - 1]) < sm) {
				if (any(dist.l > sm)) {
					id.l = tail(obs.l[dist.l > sm], 1) + 1
					x.l = knots[id.l]; x.r = knots[pos]
				} else {
					x.l = knots[1]; x.r = knots[pos]
				}
				
			}
			if (dist.r[1] < (-sm)) {
				x.r = knots[pos]
			} else if (abs(dist.r[1]) < sm) {
				if (any(dist.r < (-sm))) {
					id.r = (obs.r[dist.r < (-sm)])[1] - 1
					#x.l = knots[pos];
					x.r = knots[id.r]
				} else {
					#x.l = knots[pos]; 
					x.r = knots[m2]
				}
			}
		}
		chpt = (x.l + x.r) / 2
	} 
	if (dist[pos] < (-sm)) {
		if (pos == 1) {
			chpt = knots[pos]
		} else {
			if (any(dist > sm)) {
				if (dist[pos - 1] > sm) {
					x.l = knots[pos - 1]; x.r = knots[pos]
					y.l = dist[pos - 1]; y.r = dist[pos]
					sub = c(x.l, x.r)
				} 
				if (abs(dist[pos - 1]) < sm) {
					obs.l = 1:(pos - 1)
					dist.l = dist[obs.l]
					id.l = obs.l[dist.l > sm]
					x.l = knots[tail(id.l, 1) + 1]; x.r = knots[pos - 1]
					chpt = (x.l + x.r) / 2
					#sub = c(x.l, x.r)
					y.l = dist[tail(id.l, 1) + 1]; y.r = dist[pos - 1]
				}
			} else if (all(abs(dist[1:(pos - 1)]) < sm)) {
				x.l = knots[1]; x.r = knots[pos - 1]
				chpt = (x.l + x.r) / 2
				y.l = dist[1]; y.r = dist[pos]
				#sub = c(x.l, x.r)
			}
		}
	}
	if (!is.null(sub)) {
    		a = (y.r - y.l) / (x.r - x.l)
	  	b = y.r - a * x.r 
	  	f = function(x) a * x + b
    		chpt = uniroot(f, c(x.l, x.r), tol = 1e-8)$root
	}
	if (sh == -1) {
		thb = -thb
		thb_x = -thb_x
		bhat = -bhat
	}
	if (!is.null(zmat)) {
		zcoefs = bhat[(m2 + 2):(m2 + 2 + np - 1)]
		df_obs = sum(abs(bhat) > 0)
		prior.w = 1:n*0 + 1
		w = diag(as.vector(prior.w / deriv.fun(thb, fml = family$family)))
		sse1 = sum(prior.w * (y - thb)^2)
     		if ((n - np - 1.5 * (df_obs - np)) <= 0) {
		 	  sdhat2 = sse1
     		} else {
			 sdhat2 = sse1 / (n - np - 1.5 * (df_obs - np))
     		}
     		if (wt.iter) {
			  se2 = solve(t(zmat) %*% w %*% zmat)
     		} else {
			  se2 = solve(t(zmat) %*% diag(prior.w) %*% zmat) * sdhat2
     		}		 			
     		se.beta = sqrt(diag(se2))
     		tstat = zcoefs / se.beta
     		pvals.beta = 2 * (1 - pt(abs(tstat),  n - np - 1.5 * (df_obs - np))) 
  	} else {zcoefs = pvals.beta = se.beta = tstat = NULL}
	distp = dist[pos]
	chpt = chpt * (max(xs) - min(xs)) + min(xs)
	knots = knots * (max(xs) - min(xs)) + min(xs)
	if (!wt.iter) {thb_eta_x = thb_x; thb_eta = thb}
	ans = list(fhat = thb, fhat_x = thb_x, fhat_eta = thb_eta, fhat_eta_x = thb_eta_x, bhat = bhat, sse = minsse, pos = pos, chpt = chpt, knots = knots, dmat = dmat, tr = tr, trr = trr, tru = tr_use0, dist = dist, sub = sub, zcoefs = zcoefs, pvals.beta = pvals.beta, se.beta = se.beta, tval = tstat, bmat = bmat, phi = phi, sig = sig, aics = aics, aiclst = aiclst, es = es, lambda = lambda, lams = lambdas, edfs = edfs, gcvus = gcvus, lambdas_pen = lambdas_pen, rmat = rmat, aicmat = aicmat, sigmat = sigmat, trsmat = trsmat, trrsmat = trrsmat, psmat = psmat, fsmat = fsmat, phismat = phismat, pensmat = pensmat)
	ans 
}

########
#ip.kts#
########
ip.kts = function(x, y, zmat = NULL, m = NULL, knots = NULL, q = 3, pen = 0, pnt = FALSE, sh = 1, fir = FALSE, wt.iter = FALSE, arp = FALSE, lambdas = NULL, pen_bt = NULL, p_bt = NULL, hs = NULL, ids = NULL, family = gaussian(), gcv = FALSE) {
	linkfun = family$linkfun
	cicfamily = CicFamily(family)
	llh.fun = cicfamily$llh.fun
	etahat.fun = cicfamily$etahat.fun
	gr.fun = cicfamily$gr.fun
	wt.fun = cicfamily$wt.fun
	zvec.fun = cicfamily$zvec.fun
	muhat.fun = cicfamily$muhat.fun
	ysim.fun = cicfamily$ysim.fun
	deriv.fun = cicfamily$deriv.fun
	dev.fun = cicfamily$dev.fun 
	xs = sort(x)	
	ord = order(x)
	y = y[ord]
	xs = sort(x)	
#new:
	#x = xs 
  	mult = (max(xs) - min(xs))
	x = (xs - min(xs)) / (max(xs) - min(xs))
	n = length(x)
	xu = unique(x)
	sm = 1e-5
	if (is.null(m)) {
		m0 = 4 + round(n^(1 / 9)) 
		m1 = m0
		m = NULL
		while (is.null(m)) { 
			pts = floor((n - 2) / (m1 - 1))
			rem_pts = (n - 2) %% (m1 - 1)
			if (pts > 2) {
				m = m1
			} else if (pts == 2) {
				if (rem_pts / (m1 - 1) >= 1) {
					m = m1
				} else {
					m1 = m1 - 1	
				}
			} else {
				m1 = m1 - 1
			}
		}
	}
	m2 = m
	ans = bcspl(x, m = m2, knots = knots, pnt = pnt)
	knots = ans$knots
	bmat = ans$bmat
	secder = ans$secder
	if (!is.null(zmat)) {
		bmat = cbind(bmat, zmat)
	    	np = ncol(zmat)
	} else {np = 0}
	m = length(bmat) / n
	qv0 = crossprod(bmat)
  	qv = qv0 
  	la = 1	
	dmat = dv0 = NULL
	edfs = NULL
	gcvus = NULL
	lambda = lambdas = lambdas_pen = NULL 
	if (pnt) {
		dmat = matrix(0, nrow = (m - q), ncol = m)
#  third-order
		if (q == 3) {
			for (i in 4:m) {
				dmat[i - 3, i - 3] = 1; dmat[i - 3, i - 2] = -3; dmat[i - 3, i - 1] = 3; dmat[i - 3, i] = -1
			}
		}
# second order
		if (q == 2) {
			for (i in 3:m) {
				dmat[i - 2, i - 2] = 1; dmat[i - 2, i - 1] = -2; dmat[i - 2, i] = 1
			}
		}
# first order
		if (q == 1) {
			for (i in 2:m) {
				dmat[i - 1, i - 1] = 1; dmat[i - 1, i] = -1
			}
		}
# zero order
		if (q == 0) {
			for (i in 1:m) {
				dmat[i, i] = 1
			}
		}     
    		dv0 = crossprod(dmat)
    		if (!wt.iter) {
       			if (!arp) {
          			#if (pen == 0) {
             			#	pen = find_pen(aims = 8, Q = qv0, B = bmat, D = dv0, PNT = pnt)
          			#}
				if (pen == 0) {
					if (m < 8) {
						if (!gcv) {
# the unconstr number of columns of bmat 
							edfs = m
						} else {
							edfs = m:(m+2)
						}
						ans_pen = find_pen(aims = edfs, Q = qv0, B = bmat, D = dv0, PNT = pnt, Y = y, D0 = dmat, GCV = gcv)
	             				pen = ans_pen$lam_use
						gcvus = ans_pen$gcvus
						lambdas_pen = ans_pen$lambdas 
					} else {
						if (!gcv) {
# the unconstr number of columns of bmat 
							edfs = 8
						} else {
							edfs = 8:10
						}
						ans_pen = find_pen(aims = edfs, Q = qv0, B = bmat, D = dv0, PNT = pnt, Y = y, D0 = dmat, GCV = gcv)
						pen = ans_pen$lam_use
						gcvus = ans_pen$gcvus
						lambdas_pen = ans_pen$lambdas 
					}
		        	}
          			lambdas = pen 
       			}
       			if (arp & is.null(p_bt)) {        
	       			if (m < 10) {
            				aims = 5:m
         			}  else {
            				aims = 5:10
         			}
				edfs = aims
         			if (is.null(lambdas)) {
					ans_pen = find_pen(aims = aims, Q = qv0, B = bmat, D = dv0, PNT = pnt, Y = y, D0 = dmat, GCV = gcv)
            				lambdas = ans_pen$lam_use
         			}
       			} 
       			if (!is.null(p_bt)) {
          			lambdas = as.vector(pen_bt)
       			}
       			la = length(lambdas)
    		}
	} 
	dfmat = bmat %*% solve(qv, t(bmat))
	tru = tr_use = sum(diag(dfmat))
	dfmat = bmat %*% solve(qv, t(bmat))
	aics = minsse = nrep = pos = pos2 = NULL
	tr = sig = aics = es = rmat = phi = NULL
	aicmat = sigmat = trsmat = trrsmat = trslst = psmat = fsmat = phismat = psmat = pensmat = NULL
 	tr_use0 = tr_use = trr_use = NULL
	if (!wt.iter) {
		if (sh == -1) {
#concave-convex
			y = -y
		}
		#if (!is.null(zmat)) {
		#	bmat = cbind(bmat, zmat)
		#	qv0 = crossprod(bmat)
	    	#	np = ncol(zmat)
		#} else {np = 0}
		cv = crossprod(bmat, y)
		smat = secder 
		if (!is.null(zmat)) {
			nilmat = matrix(0, nrow = nrow(smat), ncol = ncol(zmat))
			smat = cbind(smat, nilmat)
		}
    		tdf = ncol(bmat)
    		imat = diag(tdf)
#new
    		ans_x1 = bcsplfirderi(x[1], knots)
    		ans_xn = bcsplfirderi(x[n], knots) 
    		for (ila in 1:la) {
        		if (pnt) {
           			pen = lambdas[ila]
           			qv = qv0 + pen * dv0
				dfmat = bmat %*% solve(qv, t(bmat))
				tr_use0 = c(tr_use0, sum(diag(dfmat)))
        		} else {qv = qv0}   
        		minsse = sum(y^2)
		    	for (j in 1:m2) {
            			sl = smat #-secder
				sl[j:m2, ] = -sl[j:m2, ]
			      	if (fir) {
#increase if sh == 1; decrease if sh == -1 
					row_add = matrix(1:(2 * (ncol(sl))) * 0, nrow = 2)
				       	sl = rbind(sl, row_add)
               				sl[m2 + 1, 1:(m2 + 2)] = ans_x1 
               				sl[m2 + 2, 1:(m2 + 2)] = ans_xn 
            			}
			      	qans = qprog(qv, cv, sl, 1:nrow(sl)*0, msg = FALSE)
            			bh = fitted(qans)
            			qdf = qans$df
            			fhat = bmat %*% bh
            			fhat_x = bmat[, 1:(m2 + 2), drop = FALSE] %*% bh[1:(m2 + 2), ,drop = FALSE]
            			sse = sum((y - fhat)^2)
			      	if (sse < minsse) {thb = fhat; thb_x = fhat_x; minsse = sse; pos = j; bhat = bh; tr = tr_use0; lambda = pen}   
       			}
       			pos2 = c(pos2, pos)
       			pos_use = round(mean(pos2))
    		}
   	 	if (arp) {
        		st = pos_use - 2
        		ed = pos_use + 2   
        		if (st < 1) st = 1
        		if (ed > m2) ed = m2
        		minllh = sum(y^2) #1e+3
			aiclst = fslst = siglst = pslst = phislst = penslst = trslst = trrslst = list()
			iterlst = 0
			miniter = NULL
        		for (j in st:ed) {
				iterlst = iterlst + 1
            			sl = smat 
			      	sl[j:m2, ] = -sl[j:m2, ]
            			if (fir) {
#increase if sh == 1; decrease if sh == -1 
               				row_add = matrix(1:(2 * (ncol(sl))) * 0, nrow = 2)
				       	sl = rbind(sl, row_add)
               				sl[m2 + 1, 1:(m2 + 2)] = ans_x1 
               				sl[m2 + 2, 1:(m2 + 2)] = ans_xn  
           			}
           			qans = make_arp(y, pnt = pnt, dmat = dmat, bmat = bmat, sl = sl, p_max = 2, m2 = m2, lambdas = lambdas, pen_fit = pen_bt, p_fit = p_bt, dv0 = dv0, qv0 = qv0, cv = cv, hs = hs, ids = ids)
				bh = qans$thetahat
            			theta = qans$fhat_use
				theta_x = bmat[, 1:(m2 + 2), drop = FALSE] %*% bh[1:(m2 + 2), , drop = FALSE] 
				aiclst[[iterlst]] = qans$aics	
				siglst[[iterlst]] = qans$sigs
				trslst[[iterlst]] = qans$trs
				trrslst[[iterlst]] = qans$trrs
				fslst[[iterlst]] = qans$fhats
				phislst[[iterlst]] = qans$phis
				penslst[[iterlst]] = qans$pens 
            			r_use = qans$r_use
				phi_use = qans$phi_use
            			sig_use = qans$sig_use
            			aics_use = qans$aics
            			tr_use = qans$tr_use
				trr_use = qans$trr_use
            			es_use = qans$es_use
            			lambda_use = qans$lambda_use          								
				llh = crossprod((y - theta), chol2inv(chol(r_use))) %*% (y - theta)					
				if (llh < minllh) {thb = theta; thb_x = theta_x; minllh = llh; pos = j; bhat = bh; miniter = iterlst; phi = phi_use; sig = sig_use; aics = aics_use; tr = tr_use; trr = trr_use; es = es_use; lambda = lambda_use; rmat = r_use}
        		}                                
			aicmat = aiclst[[miniter]]
			sigmat = siglst[[miniter]] 
			trsmat = trslst[[miniter]]
			trrsmat = trrslst[[miniter]] 
			fsmat = fslst[[miniter]]   
			phismat = phislst[[miniter]] 
			pensmat = penslst[[miniter]]  
			id_minaic = which(aicmat == min(aicmat))
			aics = aicmat	                 
    		}
	} else { 
#concave-convex
        	if (pnt) {
#print (pen)
			if (pen == 0) {
				if (m < 8) {
					if (!gcv) {
# the unconstr number of columns of bmat 
						edfs = m
					} else {
						edfs = m:(m+2)
					}
					ans_pen = find_pen(aims = edfs, Q = qv0, B = bmat, D = dv0, PNT = pnt, Y = y, D0 = dmat, GCV = gcv)
        				pen = ans_pen$lam_use
					gcvus = ans_pen$gcvus
					lambdas_pen = ans_pen$lambdas 
				} else {
					if (!gcv) {
# the unconstr number of columns of bmat 
						edfs = 8
					} else {
						edfs = 8:10
					}
					ans_pen = find_pen(aims = edfs, Q = qv0, B = bmat, D = dv0, PNT = pnt, Y = y, D0 = dmat, GCV = gcv)
					pen = ans_pen$lam_use
					gcvus = ans_pen$gcvus
					lambdas_pen = ans_pen$lambdas 
				}
			} 
			qv = qv0 + pen * dv0
			lambda = pen; lambdas = pen
         	} 
         	dfmat = bmat %*% solve(qv, t(bmat))
         	tr_use = sum(diag(dfmat))
         	tru = tr_use
                tdf = ncol(bmat)
                imat = diag(tdf)
         	if (family$family == "binomial") {
# use solve.QP if the last b > 1
# no zmat; no vmat
              		#cv = crossprod(bmat, y)
		        smat = secder 
              		#minsse = sum(y^2)
              		minllh = 1e+3
		        for (j in 1:m2) {
#print (j)
				sl = smat 
			        sl[j:m2, ] = -sl[j:m2, ]            
                  		#if (sh == 1) {
                     		row_add = matrix(1:(3 * ncol(sl)) * 0, nrow = 3)
                     		sl = rbind(sl, row_add)
                     		sl[m2 + 1, 1] = 1
                     		sl[m2 + 2, 1:(m2 + 2)] = bcsplfirderi(x[1], knots = knots)
                     		sl[m2 + 3, 1:(m2 + 2)] = bcsplfirderi(x[n], knots = knots)
                  		#}
                  		diff = 1
                  		nrep = 0 
                  		sm = 1e-5
                  		muhat = 1:n*0 + .5
                  		wt = 1 / muhat / (1 - muhat)
                  		while (diff > sm & nrep < 100){
#print (nrep)
                        		oldmu = muhat	
                        		nrep = nrep + 1										
                        		if (pnt) {
                           			qv = crossprod(bmat, diag(wt)) %*% bmat + pen * dv0    
                        		} else {
                          			qv = crossprod(bmat, diag(wt)) %*% bmat
                        		}
                        		cv = crossprod(bmat, diag(wt)) %*% y
                        		ans = qprog(qv, cv, sl, 1:nrow(sl)*0, msg = FALSE)
                        		bh = fitted(ans)
					if (bh[m2 + 2] > 1) {
                           			B = matrix(1:ncol(sl)*0, nrow = 1, ncol = ncol(sl))  
                           			B[m2 + 2] = 1
                           			b0 = qr.solve(qr(B), 1)
                           			H = qr.Q(qr(t(B)), complete = TRUE)[, -(1:(qr(t(B))$rank)), drop = FALSE]	
                           			atil = sl %*% H
                           			bvec = -sl %*% b0
                           			if (pnt) {
							qv1 = t(H) %*% (crossprod(bmat, diag(wt)) %*% bmat + pen * dv0) %*% H 
                           			} else {
                             				qv1 = t(H) %*% crossprod(bmat, diag(wt)) %*% bmat %*% H
                           			}
                           			cv1 = t(H) %*% crossprod(bmat, diag(wt)) %*% (y - bmat %*% b0)
#new: scale qv1 and cv1				
						const = norm(qv1, "2")
                           			qans = solve.QP(qv1 / const, cv1 / const, t(atil), bvec)
                           			phi = qans$solution
                           			btil = H %*% phi                        			
						bh = btil + b0
                        		}   
                        		muhat = bmat %*% bh
	                      		diff = mean((muhat - oldmu)^2)	
		                    	if (any(round(muhat, 8) == 0)| any(round(muhat, 8) == 1)) {
			                     wt = rep(1e+2, n)
			                     id = round(muhat, 8) > 0 & round(muhat, 8) < 1
			                     wt[id] = 1 / muhat[id] / (1 - muhat[id])
                        		} else {wt = as.vector(1 / muhat / (1 - muhat))} 
                  		}
#tr when mu_n > 1???
                     		d_use = t(ans$xmat)
                     		qdf = ans$df
                     		if (qdf < tdf) {
                       		 	pd = d_use %*% solve(crossprod(d_use), t(d_use))
                        		pmat0 = imat - pd
                     		} else {
                       			pmat0 = imat
                     		} 
                     		umat = chol(qv)
                     		iumat = diag(ncol(umat))
                     		uinv = backsolve(umat, iumat)
                     		bu = bmat %*% uinv
                     		pmat = bu %*% tcrossprod(pmat0, bu) %*% diag(wt)   
                     		tr_use = sum(diag(pmat))  
                  		muhat[round(muhat, 8) == 0] = sm
		        	muhat[round(muhat, 8) == 1] = 1
                  		llh = llh.fun(y, muhat, etahat = NULL, n, weights = NULL, fml = family$family)
                  		if (llh < minllh) {thb = muhat; thb_x = muhat; minllh = llh; pos = j; bhat = bh; tr = tr_use}
             		}               
      		}
      		if (family$family == "poisson") {
         		smat = secder 
	 		minllh = 1e+3
         		for (j in 1:m2) {
             			sl = smat 
             			sl[j:m2, ] = -sl[j:m2, ] 
             			#if (sh == 1) {
             			row_add = matrix(1:(3 * ncol(sl)) * 0, nrow = 3)
             			sl = rbind(sl, row_add)
             			sl[m2 + 1, 1] = 1
             			sl[m2 + 2, 1:(m2 + 2)] = bcsplfirderi(x[1], knots = knots)
             			if (fir) {
             				sl[m2 + 3, 1:(m2 + 2)] = bcsplfirderi(x[n], knots = knots)
             			} else {
             				sl[m2 + 3, m2 + 2] = 1
            			} 
             			#}
             			diff = 1
             			nrep = 0 
             			sm = 1e-5
             			muhat = 1:n*0 + mean(y)
             			etahat = 1:n*0
             			wt = 1 / muhat 
             			while (diff > sm & nrep < 100){
                   			oldmu = muhat	
                   			nrep = nrep + 1										
                  			 if (pnt) {
                      				qv = crossprod(bmat, diag(wt)) %*% bmat + pen * dv0    
                   			} else {
                     				qv = crossprod(bmat, diag(wt)) %*% bmat
                   			}
                   			cv = crossprod(bmat, diag(wt)) %*% y
                   			ans = qprog(qv, cv, sl, 1:nrow(sl)*0, msg = FALSE)
                   			bh = fitted(ans)
                   			muhat = bmat %*% bh
                   			diff = mean((muhat - oldmu)^2)	
                   			if (any(round(muhat, 8) == 0)) {
                     				wt = rep(1e+2, n)
                      				id = round(muhat, 8) > 0
                      				wt[id] = 1 / muhat[id]
                   			} else {wt = as.vector(1 / muhat)} 
            			 }
             			d_use = t(ans$xmat)
             			qdf = ans$df
            			if (qdf < tdf) {
               				pd = d_use %*% solve(crossprod(d_use), t(d_use))
               				pmat0 = imat - pd
             			} else {
               				pmat0 = imat
             			} 
             			umat = chol(qv)
             			iumat = diag(ncol(umat))
             			uinv = backsolve(umat, iumat)
             			bu = bmat %*% uinv
             			pmat = bu %*% tcrossprod(pmat0, bu) %*% diag(wt)   
             			tr_use = sum(diag(pmat))  
             			muhat[round(muhat, 8) == 0] = sm
             			etahat = log(muhat)
             			llh = llh.fun(y, muhat, etahat, n, weights = NULL, fml = family$family)
             			if (llh < minllh) {thb = muhat; thb_x = muhat; minllh = llh; pos = j; bhat = bh; tr = tr_use}
        		}   
        		thb[round(thb, 8) == 0] = 0                    
      		}
	}
  	dist = round(secder %*% bhat[1:(m2 + 2), ,drop = FALSE], 6) 
	sub = iter = NULL
	if (abs(dist[pos]) < sm) {
		if (pos == 1) {
			obs.r = 2:m2
			dist.r = dist[obs.r]
			if (dist.r[1] < 0) {
				x.l = knots[pos]; x.r = knots[pos]
			} else if (abs(dist.r[1]) < sm) {
				if (any(dist.r < 0)) {
					id.r = (obs.r[dist.r < 0])[1] - 1
					x.l = knots[pos]; x.r = knots[id.r]
				} else {
					x.l = knots[pos]; x.r = knots[m2]
				}
			}
		} else if (pos == m2) {
			obs.l = 1:(m2 - 1)
			dist.l = dist[obs.l]
			if (dist.l[m2 - 1] > sm) {
				x.l = knots[pos]; x.r = knots[pos]
			} else if (abs(dist.l[m2 - 1]) < sm) {
				if (any(dist.l > sm)) {
					id.l = tail(obs.l[dist.l > sm], 1) + 1
					x.l = knots[id.l]; x.r = knots[pos]
				} else {
					x.l = knots[1]; x.r = knots[pos]
				}
			}
		} else {
			obs.l = 1:(pos - 1); obs.r = (pos + 1):m2
			dist.l = dist[obs.l]; dist.r = dist[obs.r]
			if (dist.l[pos - 1] > sm) {
				x.l = knots[pos]
			} else if (abs(dist.l[pos - 1]) < sm) {
				if (any(dist.l > sm)) {
					id.l = tail(obs.l[dist.l > sm], 1) + 1
					x.l = knots[id.l]; x.r = knots[pos]
				} else {
					x.l = knots[1]; x.r = knots[pos]
				}	
			}
			if (dist.r[1] < 0) {
				x.r = knots[pos]
			} else if (abs(dist.r[1]) < sm) {
				if (any(dist.r < (-sm))) {
					id.r = (obs.r[dist.r < 0])[1] - 1
					x.l = knots[pos]; x.r = knots[id.r]
				} else {
					x.l = knots[pos]; x.r = knots[m2]
				}
			}
		}
		chpt = (x.l + x.r) / 2
	} 
	if (dist[pos] < (-sm)) {
		if (pos == 1) {
			chpt = knots[pos]
		} else {
			if (any(dist > sm)) {
				if (dist[pos - 1] > sm) {
					x.l = knots[pos - 1]; x.r = knots[pos]
					y.l = dist[pos - 1]; y.r = dist[pos]
					sub = c(x.l, x.r)
				} 
				if (abs(dist[pos - 1]) < sm) {
					obs.l = 1:(pos - 1)
					dist.l = dist[obs.l]
					id.l = obs.l[dist.l > sm]
					x.l = knots[tail(id.l, 1) + 1]; x.r = knots[pos - 1]
					chpt = (x.l + x.r) / 2
					#sub = c(x.l, x.r)
					y.l = dist[tail(id.l, 1) + 1]; y.r = dist[pos - 1]
				}
			} else if (all(abs(dist[1:(pos - 1)]) < sm)) {
				x.l = knots[1]; x.r = knots[pos - 1]
				chpt = (x.l + x.r) / 2
				y.l = dist[1]; y.r = dist[pos]
				#sub = c(x.l, x.r)
			}
		}
	}	
#print (sub)
	if (!is.null(sub)) {
		a = (y.r - y.l) / (x.r - x.l)
	    	b = y.r - a * x.r 
	    	f = function(x) a * x + b
      		chpt = uniroot(f, c(x.l, x.r), tol = 1e-8)$root
   	}
   	if (sh == -1 & !wt.iter) {
		thb = -thb
		thb_x = -thb_x
		bhat = -bhat
   	} 
	if (!is.null(zmat)) {
		zcoefs = bhat[(m2 + 3):(m2 + 3 + np - 1)]
     		df_obs = sum(abs(bhat) > 0)
		#df_obs = df		
		prior.w = 1:n*0 + 1
 		w = diag(as.vector(prior.w / deriv.fun(thb, fml = family$family)))
     		sse1 = sum(prior.w * (y - thb)^2)
     		if ((n - np - 1.5 * (df_obs - np)) <= 0) {
			sdhat2 = sse1
     		} else {
			sdhat2 = sse1 / (n - np - 1.5 * (df_obs - np))
     		}
     		if (wt.iter) {
			se2 = solve(t(zmat) %*% w %*% zmat)
     		} else {
			se2 = solve(t(zmat) %*% diag(prior.w) %*% zmat) * sdhat2
     		}		 			
     		se.beta = sqrt(diag(se2))
     		tstat = zcoefs / se.beta
     		pvals.beta = 2 * (1 - pt(abs(tstat),  n - np - 1.5 * (df_obs - np))) 
  	} else {zcoefs = pvals.beta = se.beta = tstat = NULL}
	chpt = chpt * (max(xs) - min(xs)) + min(xs)
#new: scale knots back
	knots = knots * (max(xs) - min(xs)) + min(xs)
	ans = list(fhat = thb, fhat_x = thb_x, bhat = bhat, sse = minsse, chpt = chpt, knots = knots, tr = tr, tru = tru, pos = pos, dist = dist, sub = sub, zcoefs = zcoefs, pvals.beta = pvals.beta, se.beta = se.beta, tval = tstat, bmat = bmat, phi = phi, sig = sig, aics = aics, es = es, lambda = lambda, lams = lambdas, edfs = edfs, gcvus = gcvus, lambdas_pen = lambdas_pen, rmat = rmat, aicmat = aicmat, sigmat = sigmat, trsmat = trsmat, trrsmat = trrsmat, fsmat = fsmat, phismat = phismat, pensmat = pensmat)
	 ans 
}

###########
#jp.pts   #
###########		
jp.pts = function(x, y, zmat = NULL, jpt, i = NULL, m = NULL, knots = NULL, q = 3, pen = 0, pnt = FALSE, up = TRUE, trd1 = -1, trd2 = -1, constr = TRUE, wt.iter = FALSE, arp = FALSE, lambdas = NULL, pen_bt = NULL, p_bt = NULL, hs = NULL, ids = NULL, family = gaussian(), gcv = FALSE) {	
	linkfun = family$linkfun
	cicfamily = CicFamily(family)
	llh.fun = cicfamily$llh.fun
	etahat.fun = cicfamily$etahat.fun
	gr.fun = cicfamily$gr.fun
	wt.fun = cicfamily$wt.fun
	zvec.fun = cicfamily$zvec.fun
	muhat.fun = cicfamily$muhat.fun
	ysim.fun = cicfamily$ysim.fun
	deriv.fun = cicfamily$deriv.fun
	dev.fun = cicfamily$dev.fun 	
	sm = 1e-5
	if (!up) {
		y = -y; trd1 = -trd1; trd2 = -trd2
	}
	xs = sort(x)	
	ord = order(x)
	y = y[ord]	
	n = length(x)
	#xu = unique(x)
	if (is.null(m)) {
		m0 = 4 + round(n^(1 / 9)) 
		m1 = m0
		m = NULL
		while (is.null(m)) { 
			pts = floor((n - 2) / (m1 - 1))
			rem_pts = (n - 2) %% (m1 - 1)
			if (pts > 2) {
				m = m1
			} else if (pts == 2) {
				if (rem_pts / (m1 - 1) >= 1) {
					m = m1
				} else {
					m1 = m1 - 1	
				}
			} else {
				m1 = m1 - 1
			}
		}
	}
	ans = bqspl(x, m = m, knots = knots, pic = FALSE, pnt = pnt)
	knots = ans$knots
	slopes = ans$slopes
	bmat = ans$bmat
	m = length(knots)
	kobs = 1:m
	m0 = length(bmat) / n 
	mj = m0 + 2
	obs = 1:n
	#jp is x[i]
	i = obs[xs == jpt]	
	#i = obs[round(xs, 10) == jpt]
#print (i)
#if jp is at some x[i], then there will be tiny increase
	#jp = (x[i]+x[i+1])/2
	#jp0 = x[i]
	djump = matrix(0, nrow = n, ncol = mj)
	djump[, 1:(mj - 2)] = bmat
	djump[1:i, mj - 1] = (2 * i + 1) / 2 - n; djump[(i + 1):n, mj - 1] = (2 * i + 1) / 2
	djump[1:i, mj] = 0; djump[(i + 1):n, mj] = x[x > (x[i] + x[i + 1]) / 2] - (x[i] + x[i + 1] ) / 2
	m_i = mean(djump[, mj])
	djump[, mj] = djump[, mj] - mean(djump[, mj])
#new:
	if (!is.null(zmat)) {
		djump = cbind(djump, zmat)
		np = ncol(zmat)
	}
#new: pnt, qv
	gcvus = NULL
	edfs = NULL
	lambda = lambdas = lambdas_pen = NULL
	tr = tru = sub = kts_in = NULL
	rmat = phi = sig = aics = tr = trr = es = theta = etahat = theta_x = etahat_x = NULL
	if (constr) {
		dist = knots - (x[i] + x[i + 1]) / 2
		obs = 1:m
		pos = NULL
		if (any(round(dist, 8) == 0)) {
			pos = min(obs[round(dist, 8) == 0])
			sub = c(knots[pos-1], knots[pos+1])
		} else {
			for (im in 1:(m - 1)) {
				if (dist[im] < 0 & dist[im+1] > 0) {
					pos = im
					break 
				}
			}
			sub = c(knots[pos], knots[pos+1])
		}
		bool = knots >= sub[1] & knots <= sub[2]
		kts_in = knots[bool]
		id_in = kobs[bool]
		kn = 1:(m + 3) < 0
		bool1 = knots > (x[i]+x[i+1])/2 & knots <= sub[2]		
		kn[1:m] = bool1

		smat = matrix(0, nrow = m + 3, ncol = mj)
		if (trd1 == -1 & trd2 == -1) {
			smat[id_in, 1:(mj - 2)] = -slopes[id_in, ]
#useless case:
		} else if (trd1 == 1 & trd2 == 1) {
			smat[id_in, 1:(mj - 2)] = slopes[id_in, ]
		} else if (trd1 * trd2 < 0) {
			dist = knots - (x[i] + x[i+1]) / 2
			obs = 1:m
			pos = NULL
			if (any(round(dist, 8) == 0)) {
				pos = min(obs[round(dist, 8) == 0])
			} else {
				for (im in 1:(m - 1)) {
					if (dist[im] < 0 & dist[im+1] > 0) {
						pos = im
						break 
					}
				}
			}
			bool2 = knots >= sub[1] & knots <= knots[pos]
			kts_in2 = knots[bool2]
			id_in2 = kobs[bool2]

			bool3 = knots > knots[pos] & knots <= sub[2]
			kts_in3 = knots[bool3]
			id_in3 = kobs[bool3]

			smat[id_in, 1:(mj - 2)] = slopes[id_in, ]
			if (trd1 == - 1 & trd2 == 1) {		
				smat[id_in2, 1:(mj - 2)] = -smat[id_in2, 1:(mj - 2)]		
			} 
			if (trd1 == 1 & trd2 == -1) {
				smat[id_in3, 1:(mj - 2)] = -smat[id_in3, 1:(mj - 2)]		
			}
		}
smat[m + 1, mj - 1] = 1
		if (trd1 == -1 & trd2 == -1) {
			smat[, mj] = 0
			if (any(bool1)) {
				smat[kn, mj] = -1
			} 
			smat[m + 2, mj] = -1
			ansi = sl_fun((x[i] + x[i+1]) / 2, knots, slopes)
			smat[m + 2, 1:m0] = -ansi 
			smat[m + 3, 1:m0] = -ansi  
		} else if (trd1 == 1 & trd2 == 1) {
			smat[, mj] = 0
			if (any(bool1)) {
				smat[kn, mj] = 1
			} 
			smat[m + 2, mj] = 1
			ansi = sl_fun((x[i] + x[i+1]) / 2, knots, slopes)
			smat[m + 2, 1:m0] = ansi 
			smat[m + 3, 1:m0] = ansi 
		} else if (trd1 == -1 & trd2 == 1) {
			smat[, mj] = 0
			if (any(bool1)) {
				smat[kn, mj] = 1
			} 
			smat[m + 2, mj] = 1
			ansi = sl_fun((x[i] + x[i+1]) / 2, knots, slopes)
			smat[m + 2, 1:m0] = ansi 
			smat[m + 3, 1:m0] = -ansi 
		} else if (trd1 == 1 & trd2 == -1) {
#smat[m + 1, mj - 1] = 1
			smat[, mj] = 0
			if (any(bool1)) {
				smat[kn, mj] = -1
			} 
			smat[m + 2, mj] = -1
			ansi = sl_fun((x[i] + x[i+1]) / 2, knots, slopes)
			smat[m + 2, 1:m0] = -ansi 
			smat[m + 3, 1:m0] = ansi 
		} 
		use = 1:(m + 3) > 0
		kp = min(knots[knots > (x[i]+x[i+1])/2])
		if (sum(x > (x[i]+x[i+1])/2 & x <= kp) == 0){use[m + 2] = FALSE}
		kp = max(knots[knots < (x[i]+x[i+1])/2])
		if (sum(x < (x[i]+x[i+1])/2 & x >= kp) == 0){use[m + 3] = FALSE}
		smat = smat[use, ]
		if (i == 1 | i == (n - 1)) {
			smat = smat[, -mj, drop = FALSE]
			djump = djump[, -mj, drop = FALSE] 
		}	
		if (!is.null(zmat)) {
			nilmat = matrix(0, nrow = nrow(smat), ncol = ncol(zmat))
			smat = cbind(smat, nilmat)
		}
		qv0 = crossprod(djump)
		qv = qv0		
		if (pnt) {
			mj2 = ncol(djump)
			dmat = matrix(0, nrow = (mj2 - q), ncol = mj2)
#  third-order
			if (q == 3) {
				for (idm in 4:mj2) {
					dmat[idm - 3, idm - 3] = 1; dmat[idm - 3, idm - 2] = -3; dmat[idm - 3, idm - 1] = 3; dmat[idm - 3, idm] = -1
				}
			}
# second order
			if (q == 2) {
				for (idm in 3:mj2) {
					dmat[idm - 2, idm - 2] = 1; dmat[idm - 2, idm - 1] = -2; dmat[idm - 2, idm] = 1
				}
			}
# first order
			if (q == 1) {
				for (idm in 2:mj2) {
					dmat[idm - 1, idm - 1] = 1; dmat[idm - 1, idm] = -1
				}
			}
# zero order
			if (q == 0) {
				for (idm in 1:mj2) {
					dmat[idm, idm] = 1
				}
			}    
    			dv0 = crossprod(dmat)
       			if (!arp) {
          			#if (pen == 0) {
             			#	pen = find_pen(aims = 6, Q = qv0, B = djump, D = dv0, PNT = pnt)
          			#}
				if (pen == 0) {
					if (mj2 < 8) {
						if (!gcv) {
# the unconstr number of columns of bmat 
							edfs = mj2
						} else {
							edfs = mj2:(mj2+2)
						}
						ans_pen = find_pen(aims = edfs, Q = qv0, B = djump, D = dv0, PNT = pnt, Y = y, D0 = dmat, GCV = gcv)
	             				pen = ans_pen$lam_use
						gcvus = ans_pen$gcvus
						lambdas_pen = ans_pen$lambdas 
					} else {
						if (!gcv) {
# the unconstr number of columns of bmat 
							edfs = 8
						} else {
							edfs = 8:10
						}
						ans_pen = find_pen(aims = edfs, Q = qv0, B = djump, D = dv0, PNT = pnt, Y = y, D0 = dmat, GCV = gcv)
						pen = ans_pen$lam_use
						gcvus = ans_pen$gcvus
						lambdas_pen = ans_pen$lambdas 
					}
		        	}
          			lambdas = pen
       			}
       			if (arp & is.null(p_bt)) {        
	        		if (mj2 < 10) {
					aims = 5:mj2
         			} else {
            				aims = 5:10
         			}
				edfs = aims
         			if (is.null(lambdas)) {
            				if (pen == 0) {
						ans_pen = find_pen(aims = aims, Q = qv0, B = djump, D = dv0, PNT = pnt, Y = y, D0 = dmat, GCV = gcv)
               					lambdas = ans_pen$lam_use
            				} else {
               					lambdas = pen
            				}
         			}
       			} 
      			if (!is.null(p_bt)) {
          			lambdas = as.vector(pen_bt)
       			}
       			la = length(lambdas)
			qv = qv + pen * dv0
		} 
		if (!wt.iter) {
			cv = crossprod(djump, y)
			if (!arp) {
	         		dfmat = djump %*% solve(qv, t(djump))
	         		tr_use = sum(diag(dfmat))
	         		tru = tr_use
				fit = qprog(qv, cv, smat, 1:nrow(smat)*0, msg = FALSE)				
				bhat = fitted(fit)		
				theta = djump %*% bhat
				nd = ncol(djump)
				if (!is.null(zmat)) {
					nd = nd - ncol(zmat)
				}
				theta_x = djump[, 1:nd, drop = FALSE] %*% bhat[1:nd, ,drop = FALSE]
				df = fit$df
			} else {
				qans = make_arp(y, pnt = pnt, dmat = dmat, bmat = djump, sl = smat, p_max = 4, m2 = m, lambdas = lambdas, pen_fit = pen_bt, p_fit = p_bt, dv0 = dv0, qv0 = qv0, cv = cv, hs = hs, ids = ids)
				bhat = qans$thetahat
	            		theta = qans$fhat_use
				nd = ncol(djump)
				if (!is.null(zmat)) {
					nd = nd - ncol(zmat)
				}
				theta_x = djump[, 1:nd, drop = FALSE] %*% bhat[1:nd, , drop = FALSE]
	            		rmat = qans$r_use
				phi = qans$phi_use
	            		sig = qans$sig_use
#print (sig)
	            		aics = qans$aics
	            		tr = qans$tr_use
				trr = qans$trr_use
	            		es = qans$es_use
	            		lambda = qans$lambda_use
			}
		} else {
			etahat = etahat.fun(n, y, fml = family$family)
			gr = gr.fun(y, etahat, weights = NULL, fml = family$family)  
			wt = wt.fun(etahat, n, weights = NULL, fml = family$family)
			cvec = wt * etahat - gr
			qvk = crossprod(djump, diag(wt)) %*% djump
	        	if (pnt) {
				#if (pen == 0) {
	              		#	pen = find_pen(aims = 6, B = djump, D = dmat, PNT = pnt)
				#}
				if (pen == 0) {
					if (mj2 < 8) {
						if (!gcv) {
# the unconstr number of columns of bmat 
							edfs = mj2
						} else {
							edfs = mj2:(mj2+2)
						}
						ans_pen = find_pen(aims = edfs, Q = qv0, B = djump, D = dv0, PNT = pnt, Y = y, D0 = dmat, GCV = gcv)
	             				pen = ans_pen$lam_use
						gcvus = ans_pen$gcvus
						lambdas_pen = ans_pen$lambdas 
					} else {
						if (!gcv) {
# the unconstr number of columns of bmat 
							edfs = 8
						} else {
							edfs = 8:10
						}
						ans_pen = find_pen(aims = edfs, Q = qv0, B = djump, D = dv0, PNT = pnt, Y = y, D0 = dmat, GCV = gcv)
						pen = ans_pen$lam_use
						gcvus = ans_pen$gcvus
						lambdas_pen = ans_pen$lambdas 
					}
		        	}
				qvk = qvk + pen * dv0
				lambda = pen; lambdas = pen        	
			}
			cveck = crossprod(djump, cvec)			
			qans = qprog(qvk, cveck, smat, 1:nrow(smat)*0, msg = FALSE)
			bhat = qans$thetahat
			etahat = djump %*% bhat
			muhat = muhat.fun(etahat, fml = family$family)
			diff = 1
			nrep = 0
			while (diff > sm & nrep < 100) {
				oldmu = muhat	
				nrep = nrep + 1
				gr = gr.fun(y, etahat, weights = NULL, fml = family$family)	
				wt = wt.fun(etahat, n, weights = NULL, fml = family$family) 
				cvec = wt * etahat - gr
				qvk = crossprod(djump, diag(wt)) %*% djump
				if (pnt) {
					qvk = qvk + pen * dv0
				}
				cveck = crossprod(djump, cvec)	
				qans = qprog(qvk, cveck, smat, 1:nrow(smat)*0, msg = FALSE)
				bhat = qans$thetahat
				etahat = djump %*% bhat
				nd = ncol(djump)
				if (!is.null(zmat)) {
					nd = nd - ncol(zmat)
				}
				muhat = muhat.fun(etahat, fml = family$family)
				diff = mean((muhat - oldmu)^2)	
#lines(x, muhat, col = nrep + 1, lty = 2)
	       		} 
	       		etahat_x = djump[, 1:nd, drop = FALSE] %*% bhat[1:nd, ,drop = FALSE]
	       		muhat_x = muhat.fun(etahat_x, fml = family$family)
			theta = muhat
			theta_x = muhat_x
		}
	} else {
		if (!wt.iter) {
			if (i == 1 | i == (n - 1)) {
				djump = djump[, -mj, drop = FALSE] 
			}	
			qv0 = crossprod(djump)
			qv = qv0
			cv = crossprod(djump, y)
			if (pnt) {
				mj2 = ncol(djump)
				dmat = matrix(0, nrow = (mj2 - q), ncol = mj2)
#  third-order
				if (q == 3) {
					for (idm in 4:mj2) {
						dmat[idm - 3, idm - 3] = 1; dmat[idm - 3, idm - 2] = -3; dmat[idm - 3, idm - 1] = 3; dmat[idm - 3, idm] = -1
					}
				}
# second order
				if (q == 2) {
					for (idm in 3:mj2) {
						dmat[idm - 2, idm - 2] = 1; dmat[idm - 2, idm - 1] = -2; dmat[idm - 2, idm] = 1
					}
				}
# first order
				if (q == 1) {
					for (idm in 2:mj2) {
						dmat[idm - 1, idm - 1] = 1; dmat[idm - 1, idm] = -1
					}
				}
# zero order
				if (q == 0) {
					for (idm in 1:mj2) {
						dmat[idm, idm] = 1
					}
				}    
    				dv0 = crossprod(dmat)
       				if (!arp) {
					if (pen == 0) {
						if (mj2 < 8) {
							if (!gcv) {
# the unconstr number of columns of bmat 
								edfs = mj2
							} else {
								edfs = mj2:(mj2+2)
							}
							ans_pen = find_pen(aims = edfs, Q = qv0, B = djump, D = dv0, PNT = pnt, Y = y, D0 = dmat, GCV = gcv)
	             					pen = ans_pen$lam_use
							gcvus = ans_pen$gcvus
							lambdas_pen = ans_pen$lambdas 
						} else {
							if (!gcv) {
# the unconstr number of columns of bmat 
								edfs = 8
							} else {
								edfs = 8:10
							}
							ans_pen = find_pen(aims = edfs, Q = qv0, B = djump, D = dv0, PNT = pnt, Y = y, D0 = dmat, GCV = gcv)
							pen = ans_pen$lam_use
							gcvus = ans_pen$gcvus
							lambdas_pen = ans_pen$lambdas 
						}
		        		}
          				lambdas = pen
       				}
       				if (arp & is.null(p_bt)) {        
	        			if (mj2 < 10) {
						aims = 5:mj2
         				} else {
            					aims = 5:10
         				}
					edfs = aims
         				if (is.null(lambdas)) {
            					if (pen == 0) {
							ans_pen = find_pen(aims = aims, Q = qv0, B = djump, D = dv0, PNT = pnt, Y = y, D0 = dmat, GCV = gcv)
               						lambdas = ans_pen$lam_use
            					} else {
               						lambdas = pen
            					}
         				}
       				} 
      				if (!is.null(p_bt)) {
          				lambdas = as.vector(pen_bt)
       				}
       				la = length(lambdas)
				qv = qv + pen * dv0
			} 	
			if (!arp) {
				bhat = solve(qv, cv)				
				#bhat = solve(crossprod(djump), crossprod(djump, y))
				theta = djump %*% bhat		
				#nd = ncol(djump)
				#if (!is.null(zmat)) {
				#	nd = nd - ncol(zmat)
				#}
				#theta_x = djump[, 1:nd, drop = FALSE] %*% bhat[1:nd, ,drop = FALSE]
				#df = n - ncol(djump)
			} else {
				qans = make_arp(y, pnt = pnt, dmat = dmat, bmat = djump, sl = NULL, p_max = 4, m2 = m, lambdas = lambdas, pen_fit = pen_bt, p_fit = p_bt, dv0 = dv0, qv0 = qv0, cv = cv, hs = hs, ids = ids, constr = FALSE)
				bhat = qans$thetahat
	            		theta = qans$fhat_use
				#nd = ncol(djump)
				#if (!is.null(zmat)) {
				#	nd = nd - ncol(zmat)
				#}
				#theta_x = djump[, 1:nd, drop = FALSE] %*% bhat[1:nd, , drop = FALSE]
	            		rmat = qans$r_use
				phi = qans$phi_use
	            		sig = qans$sig_use
	            		aics = qans$aics
	            		tr = qans$tr_use
				trr = qans$trr_use
	            		es = qans$es_use
	            		lambda = qans$lambda_use
			}
			nd = ncol(djump)
			if (!is.null(zmat)) {
				nd = nd - ncol(zmat)
			}
			theta_x = djump[, 1:nd, drop = FALSE] %*% bhat[1:nd, ,drop = FALSE]
			df = n - ncol(djump)
		} else {
			etahat = etahat.fun(n, y, fml = family$family)
			gr = gr.fun(y, etahat, weights = NULL, fml = family$family)  
			wt = wt.fun(etahat, n, weights = NULL, fml = family$family)
			cvec = wt * etahat - gr
			qvk = crossprod(djump, diag(wt)) %*% djump
        		if (pnt) {
				#if (pen == 0) {
       		 	      	#	pen = find_pen(aims = 6, B = djump, D = dmat, PNT = pnt)
				#}
				if (pen == 0) {
					if (mj2 < 8) {
						if (!gcv) {
# the unconstr number of columns of bmat 
							edfs = mj2
						} else {
							edfs = mj2:(mj2+2)
						} 
						ans_pen = find_pen(aims = edfs, Q = qv0, B = djump, D = dv0, PNT = pnt, Y = y, D0 = dmat, GCV = gcv)
             					pen = ans_pen$lam_use
						gcvus = ans_pen$gcvus
						lambdas_pen = ans_pen$lambdas 
					} else {
						if (!gcv) {
# the unconstr number of columns of bmat 
							edfs = 8
						} else {
							edfs = 8:10
						}
						ans_pen = find_pen(aims = edfs, Q = qv0, B = djump, D = dv0, PNT = pnt, Y = y, D0 = dmat, GCV = gcv)
						pen = ans_pen$lam_use
						gcvus = ans_pen$gcvus
						lambdas_pen = ans_pen$lambdas 
					}
	        		}
				qvk = qvk + pen * dv0
				lambda = pen; lambdas = pen
       	 		}
			cveck = crossprod(djump, cvec)			
			bhat = solve(qvk, cveck)
			etahat = djump %*% bhat
			muhat = muhat.fun(etahat, fml = family$family)
			diff = 1
			nrep = 0
			while (diff > sm & nrep < 100) {
				oldmu = muhat	
				nrep = nrep + 1
				gr = gr.fun(y, etahat, weights = NULL, fml = family$family)	
				wt = wt.fun(etahat, n, weights = NULL, fml = family$family) 
				cvec = wt * etahat - gr
				qvk = crossprod(djump, diag(wt)) %*% djump
				if (pnt) {
					qvk = qvk + pen * dv0
				}
				cveck = crossprod(djump, cvec)	
				bhat = solve(qvk, cveck)
				etahat = djump %*% bhat
				nd = ncol(djump)
				if (!is.null(zmat)) {
					nd = nd - ncol(zmat)
				}
				muhat = muhat.fun(etahat, fml = family$family)
				diff = mean((muhat - oldmu)^2)	
       			} 
			etahat_x = djump[, 1:nd, drop = FALSE] %*% bhat[1:nd, ,drop = FALSE]
			muhat_x = muhat.fun(etahat_x, fml = family$family)
			theta = muhat
			theta_x = muhat_x
		}
	}
	if (!up) {
		theta = -theta
		theta_x = -theta_x
		bhat = -bhat
	}
	if (!is.null(zmat)) {
		mj2 = ncol(djump)
		zcoefs = bhat[(mj2 - np + 1):mj2]
     		df_obs = sum(abs(bhat) > 0)
		#df_obs = df		
		prior.w = 1:n*0 + 1
 		w = diag(as.vector(prior.w / deriv.fun(theta, fml = family$family)))
     		sse1 = sum(prior.w * (y - theta)^2)
     		if ((n - np - 1.5 * (df_obs - np)) <= 0) {
			sdhat2 = sse1
     		} else {
			sdhat2 = sse1 / (n - np - 1.5 * (df_obs - np))
     		}
     		if (wt.iter) {
			se2 = solve(t(zmat) %*% w %*% zmat)
     		} else {
			se2 = solve(t(zmat) %*% diag(prior.w) %*% zmat) * sdhat2
     		}		 			
     		se.beta = sqrt(diag(se2))
     		tstat = zcoefs / se.beta
     		pvals.beta = 2 * (1 - pt(abs(tstat),  n - np - 1.5 * (df_obs - np))) 
  	} else {zcoefs = pvals.beta = se.beta = tstat = NULL}
	sse = sum((y - theta)^2)
	ans = list(pos = i, chpt = jpt, fhat = theta, fhat_x = theta_x, fhat_eta = etahat, fhat_eta_x = etahat_x, df = df, sse = sse, knots = knots, sub = sub, kts_in = kts_in, bhat = bhat, zcoefs = zcoefs, m_i = m_i, sub = sub, lams = lambdas, lambda = lambda, edfs = edfs, gcvus = gcvus, lambdas_pen = lambdas_pen, tr = tr, trr = trr, tru = tru, phi = phi, sig = sig, aics = aics, es = es, lambda = lambda, lams = lambdas, rmat = rmat)
	ans
}

############################################################
#make the first derivative at a jump point by interpolation#
############################################################
sl_fun = function(xinterp, knots, slopes) {
	nk = length(knots)
	obs = 1:nk
	dist = knots - xinterp
	id = NULL; id1 = NULL
	if (any(round(dist, 8) == 0)) {
		id = min(obs[round(dist, 8) == 0])
	} else {	
		for (i in 1:(nk - 1)) {
			if (dist[i] < 0 && dist[i + 1] > 0) {
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

#######################
#bqspline basis matrix# 
#######################
bqspl = function(x, m = NULL, knots = NULL, pic = FALSE, spl = TRUE, x0 = NULL, pnt = FALSE) {
############
#make knots#  
############
	bmat = slopes = NULL
	xu = unique(x)
#knots are equal x quantile of length = m spaced on the support of min(x) to max(x)
  	if (!is.null(knots)) {
		m = length(knots)
	} else {
		if (!is.null(m)) { 
			#knots = quantile(xu, probs = seq(0, 1, length = m), names = FALSE)
			knots = 0:(m-1)/(m-1)*(max(x)-min(x))+min(x)
			#if (pnt) {
			#	knots = 0:(m-1)/(m-1)*(max(xu)-min(xu))+min(xu)
			#} else {
			#	knots = quantile(xu, probs = seq(0, 1, length = m), names = FALSE)
			#}
		} else {
			m0 = 4 + round(n^(1 / 9)) 
			m1 = m0
			#m = NULL
			while (is.null(m)) { 
				pts = floor((n - 2) / (m1 - 1))
				rem_pts = (n - 2) %% (m1 - 1)
				if (pts > 3) {
					m = m1
				} else if (pts == 3) {
					if (rem_pts / (m1 - 1) >= 1) {
						m = m1
					} else {
						m1 = m1 - 1	
					}
				} else {
					m1 = m1 - 1
				}
			}
			#knots = quantile(xu, probs = seq(0, 1, length = m), names = FALSE)
		 	knots = 0:(m-1)/(m-1)*(max(x)-min(x))+min(x)
			#if (pnt) {
			#	knots = 0:(m-1)/(m-1)*(max(xu)-min(xu))+min(xu)
			#} else {
			#	knots = quantile(xu, probs = seq(0, 1, length = m), names = FALSE)
			#}
		}
	}
	t = 1:(m+4)*0
#new: used for predict
if (is.null(x0)) {
	t[1] = t[2] = min(x)
	t[m+3] = t[m+4] = max(x)
} else {
	t[1] = t[2] = min(x0)
	t[m+3] = t[m+4] = max(x0)
}
	t[3:(m+2)] = knots
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


#######################
#bcspline basis matrix#
#######################
bcspl = function(x, m = NULL, knots = NULL, pic = FALSE, pnt = FALSE, spl = TRUE, x0 = NULL) {
############
#make knots#  
############
	bmat = secder = NULL
	n = length(x) #?length(xu)
	xu = unique(x)	
	if (!is.null(knots)) {
		m = length(knots)
	} else {
		if (!is.null(m)) { 
			knots = 0:(m-1)/(m-1)*(max(x)-min(x))+min(x)
			#knots = quantile(xu, probs = seq(0, 1, length = m), names = FALSE)
			#if (pnt) {
			#	knots = 0:(m-1)/(m-1)*(max(xu)-min(xu))+min(xu)
			#} else {
			#	knots = quantile(xu, probs = seq(0, 1, length = m), names = FALSE)
			#}
		} else {
			m0 = 4 + round(n^(1 / 9)) 
			m1 = m0
			#m = NULL
			while (is.null(m)) { 
				pts = floor((n - 2) / (m1 - 1))
				rem_pts = (n - 2) %% (m1 - 1)
				if (pts > 3) {
					m = m1
				} else if (pts == 3) {
					if (rem_pts / (m1 - 1) >= 1) {
						m = m1
					} else {
						m1 = m1 - 1	
					}
				} else {
					m1 = m1 - 1
				}
			}
			knots = 0:(m-1)/(m-1)*(max(x)-min(x))+min(x)
			#knots = quantile(xu, probs = seq(0, 1, length = m), names = FALSE)
			#if (pnt) {
			#	knots = 0:(m-1)/(m-1)*(max(xu)-min(xu))+min(xu)
			#} else {
			#	knots = quantile(xu, probs = seq(0, 1, length = m), names = FALSE)
			#}
		}
	}
	t = 1:(m+6)*0
#new: used for predict
if (is.null(x0)) {
	t[1:3] = min(x)
	t[(m+4):(m+6)] = max(x)
} else {
  t[1:3] = min(x0)
	t[(m+4):(m+6)] = max(x0)
}
	t[4:(m+3)] = knots
##########################
#make cubic bspline basis#
##########################
	#n = length(x)
if (spl) {
	bmat = matrix(0, nrow = n, ncol = m+2)
#1st edge
	i = 1; j = i+1
	bool = x < t[5] & x >= t[4]
	bmat[bool, 1] = (t[5] - x[bool]) / (t[5] - t[2]) * (t[j+3] - x[bool])**2 / (t[j+3] - t[j+1]) / (t[j+3] - t[j+2])
	
#2nd edge 1st part
	i = 2; j = i+1
	bool = x <= t[5] & x >= t[4]
	a = (t[i+3] - x[bool])**2 / (t[i+3] - t[i+1]) / (t[i+3] - t[i+2])
	b = (x[bool] - t[j]) * (t[j+2] - x[bool]) / (t[j+2] - t[j]) / (t[j+2] - t[j+1]) + (t[j+3] - x[bool]) * (x[bool] - t[j+1]) / (t[j+3] - t[j+1]) / (t[j+2] - t[j+1])
	bmat[bool, 2] = (x[bool] - t[i]) / (t[i+3] - t[i]) * a + (t[i+4] - x[bool]) / (t[i+4] - t[i+1]) * b

#2nd edge 2nd part
	bool = x <= t[6] & x > t[5]
	b = (t[j+3] - x[bool])**2 / (t[j+3] - t[j+1]) / (t[j+3] - t[j+2])
	bmat[bool, 2] = (t[i+4] - x[bool]) / (t[i+4] - t[i+1]) * b

#3rd edge 1st part
	i = 3; j = i+1
	bool = x <= t[5] & x >= t[4]
	a = (x[bool] - t[i]) * (t[i+2] - x[bool]) / (t[i+2] - t[i]) / (t[i+2] - t[i+1]) + (t[i+3] - x[bool]) * (x[bool] - t[i+1]) / (t[i+3] - t[i+1]) / (t[i+2] - t[i+1])
	b = (x[bool] - t[j])**2 / (t[j+2] - t[j]) / (t[j+1] - t[j])
	bmat[bool, 3] = (x[bool] - t[i]) / (t[i+3] - t[i]) * a + (t[i+4] - x[bool]) / (t[i+4] - t[i+1]) * b

#3rd edge 2nd part
	bool = x <= t[6] & x > t[5]
	a = (t[i+3] - x[bool])**2 / (t[i+3] - t[i+1]) / (t[i+3] - t[i+2])
	b =  (x[bool] - t[j]) * (t[j+2] - x[bool]) / (t[j+2] - t[j]) / (t[j+2] - t[j+1]) + (t[j+3] - x[bool]) * (x[bool] - t[j+1]) / (t[j+3] - t[j+1]) / (t[j+2] - t[j+1])
	bmat[bool, 3] = (x[bool] - t[i]) / (t[i+3] - t[i]) * a + (t[i+4] - x[bool]) / (t[i+4] - t[i+1]) * b

#3rd edge 3rd part
	bool = x <= t[7] & x > t[6]
	b = (t[j+3] - x[bool])**2 / (t[j+3] - t[j+1]) / (t[j+3] - t[j+2])
	bmat[bool, 3] = (t[i+4] - x[bool]) / (t[i+4] - t[i+1]) * b

#4th to (m-1)th edge
	for (i in 4:(m-1)) {
#1st part
		j = i+1
		bool = x <= t[i+1] & x >= t[i]
		a =  (x[bool] - t[i])**2 / (t[i+2] - t[i]) / (t[i+1] - t[i])
		bmat[bool, i] = (x[bool] - t[i]) / (t[i+3] - t[i]) * a

#2nd part
		bool = x <= t[i+2] & x > t[i+1]
		a = (x[bool] - t[i]) * (t[i+2] - x[bool]) / (t[i+2] - t[i]) / (t[i+2] - t[i+1]) + (t[i+3] - x[bool]) * (x[bool] - t[i+1]) / (t[i+3] - t[i+1]) / (t[i+2] - t[i+1])
		b = (x[bool] - t[j])**2 / (t[j+2] - t[j]) / (t[j+1] - t[j])
		bmat[bool, i] = (x[bool] - t[i]) / (t[i+3] - t[i]) * a + (t[i+4] - x[bool]) / (t[i+4] - t[i+1]) * b

#3rd part
		bool = x <= t[i+3] & x > t[i+2]
		a =  (t[i+3] - x[bool])**2 / (t[i+3] - t[i+1]) / (t[i+3] - t[i+2])
		b = (x[bool] - t[j]) * (t[j+2] - x[bool]) / (t[j+2] - t[j]) / (t[j+2] - t[j+1]) + (t[j+3] - x[bool]) * (x[bool] - t[j+1]) / (t[j+3] - t[j+1]) / (t[j+2] - t[j+1])
		bmat[bool, i] = (x[bool] - t[i]) / (t[i+3] - t[i]) * a + (t[i+4] - x[bool]) / (t[i+4] - t[i+1]) * b

#4th part
		bool = x <= t[i+4] & x > t[i+3]
		b =  (t[j+3] - x[bool])**2 / (t[j+3] - t[j+1]) / (t[j+3] - t[j+2])
		bmat[bool, i] = (t[i+4] - x[bool]) / (t[i+4] - t[i+1]) * b
}

#mth edge 1st part
	i = m; j = i+1
 	bool = x <= t[m+1] & x >= t[m]
	a = (x[bool] - t[i])**2 / (t[i+2] - t[i]) / (t[i+1] - t[i])
	bmat[bool, i] = (x[bool] - t[i]) / (t[i+3] - t[i]) * a

#mth edge 2nd part
	bool = x <= t[m+2] & x > t[m+1]
	a = (x[bool] - t[i]) * (t[i+2] - x[bool]) / (t[i+2] - t[i]) / (t[i+2] - t[i+1]) + (t[i+3] - x[bool]) * (x[bool] - t[i+1]) / (t[i+3] - t[i+1]) / (t[i+2] - t[i+1])
	b = (x[bool] - t[j])**2 / (t[j+2] - t[j]) / (t[j+1] - t[j])
	bmat[bool, i] = (x[bool] - t[i]) / (t[i+3] - t[i]) * a + (t[i+4] - x[bool]) / (t[i+4] - t[i+1]) * b

#mth edge 3rd part
	bool = x <= t[m+3] & x > t[m+2]
	a = (t[i+3] - x[bool])**2 / (t[i+3] - t[i+1]) / (t[i+3] - t[i+2])
	b = (x[bool] - t[j]) * (t[j+2] - x[bool]) / (t[j+2] - t[j]) / (t[j+2] - t[j+1]) + (t[j+3] - x[bool]) * (x[bool] - t[j+1]) / (t[j+3] - t[j+1]) / (t[j+2] - t[j+1])
	bmat[bool, i] = (x[bool] - t[i]) / (t[i+3] - t[i]) * a + (t[i+4] - x[bool]) / (t[i+4] - t[i+1]) * b

#(m+1)th edge 1st part
	i = m+1; j = i+1
	bool = x <= t[m+2] & x >= t[m+1]
	a = (x[bool] - t[i])**2 / (t[i+2] - t[i]) / (t[i+1] - t[i])
	bmat[bool, i] = (x[bool] - t[i]) / (t[i+3] - t[i]) * a 

#(m+1)th edge 2nd part
	bool = x <= t[m+3] & x > t[m+2]
	a = (x[bool] - t[i]) * (t[i+2] - x[bool]) / (t[i+2] - t[i]) / (t[i+2] - t[i+1]) + (t[i+3] - x[bool]) * (x[bool] - t[i+1]) / (t[i+3] - t[i+1]) / (t[i+2] - t[i+1])
	b = (x[bool] - t[j])**2 / (t[j+2] - t[j]) / (t[j+1] - t[j])
	bmat[bool, i] = (x[bool] - t[i]) / (t[i+3] - t[i]) * a + (t[i+4] - x[bool]) / (t[i+4] - t[i+1]) * b

#(m+2)th edge
	i = m+2
	bool = x <= t[m+3] & x >= t[m+2]
	a = (x[bool] - t[i])**2 / (t[i+2] - t[i]) / (t[i+1] - t[i])
	bmat[bool, m+2] = (x[bool] - t[i]) / (t[i+3] - t[i]) * a

#####################################
#make bcsplines' 2nd derivatives    #
#at knots. row: knots; column: basis# 
#####################################
	secder = matrix(0, nrow = m, ncol = (m+2))

#1st edge: t4 and t5
	secder[1,1] = 6 / (t[5] - t[2]) / (t[5] - t[3])
	#secder[2,1] = 0 

#2nd edge: t4, t5 and t6
	secder[1,2] = -6 * (1 / (t[5] - t[2]) + 1 / (t[6] - t[3])) / (t[5] - t[3])
	secder[2,2] = 6 / (t[6] - t[3]) / (t[6] - t[4])
	#secder[3,2] = 0 

#3rd edge: t4, t5, t6 and t7
	secder[1,3] = 6 / (t[6] - t[3]) / (t[5] - t[3])	
	secder[2,3] = -6 * (1 / (t[6] - t[3]) + 1 / (t[7] - t[4])) / (t[6] - t[4])
	secder[3,3] = 6 / (t[7] - t[4]) / (t[7] - t[5])
	#secder[4,3] = 0

#4th  to (m-1)th edge: t[i], t[i+1], t[i+2], t[i+3] and t[i+4]
if (m > 4) {
	for (i in 4:(m-1)) {
		#secder[i-3,i] = 0
		secder[i-2,i] = 6 / (t[i+3] - t[i]) / (t[i+2] - t[i])
		secder[i-1,i] = -6 * (1 / (t[i+3] - t[i]) + 1 / (t[i+4] - t[i+1])) / (t[i+3] - t[i+1])
		secder[i,i] = 6 / (t[i+4] - t[i+1]) / (t[i+4] - t[i+2])
		#secder[i+1,i] = 0
	}
}
#mth edge: t[m], t[m+1], t[m+2] and t[m+3]
	#secder[m-3,m] = 0
	secder[m-2,m] = 6 / (t[m+3] - t[m]) / (t[m+2] - t[m])
	secder[m-1,m] = -6 * (1 / (t[m+3] - t[m]) + 1 / (t[m+4] - t[m+1])) / (t[m+3] - t[m+1])
	secder[m,m] = 6 / (t[m+4] - t[m+2]) / (t[m+4] - t[m+1])

#(m+1)th edge: t[m+1], t[m+2] and t[m+3]
	#secder[m-2,m+1] = 0
	secder[m-1,m+1] = 6 / (t[m+4] - t[m+1]) / (t[m+3] - t[m+1])
	secder[m,m+1] = -6 * (1 / (t[m+4] - t[m+1]) + 1 / (t[m+5] - t[m+2])) / (t[m+4] - t[m+2])

#(m+2)th edge: t[m+2] and t[m+3]
	#secder[m-1,m+2] = 0
	secder[m,m+2] = 6 / (t[m+5] - t[m+2]) / (t[m+4] - t[m+2]) 
}
########################
##plot the spline basis# 
########################
xpl = bpl = NULL
if (pic == TRUE) {
	xpl = 0:1000/1000*(max(x)-min(x))+min(x)
	bpl = matrix(0, nrow = 1001, ncol = m+2)
#1st edge
	i = 1; j = i+1
	bool = xpl < t[5] & xpl >= t[4]
	bpl[bool, 1] = (t[5] - xpl[bool]) / (t[5] - t[2]) * (t[j+3] - xpl[bool])**2 / (t[j+3] - t[j+1]) / (t[j+3] - t[j+2])

#2nd edge 1st part
	i = 2; j = i+1
	bool = xpl <= t[5] & xpl > t[4]
	a = (t[i+3] - xpl[bool])**2 / (t[i+3] - t[i+1]) / (t[i+3] - t[i+2])
	b = (xpl[bool] - t[j]) * (t[j+2] - xpl[bool]) / (t[j+2] - t[j]) / (t[j+2] - t[j+1]) + (t[j+3] - xpl[bool]) * (xpl[bool] - t[j+1]) / (t[j+3] - t[j+1]) / (t[j+2] - t[j+1])
	bpl[bool, 2] = (xpl[bool] - t[i]) / (t[i+3] - t[i]) * a + (t[i+4] - xpl[bool]) / (t[i+4] - t[i+1]) * b

#2nd edge 2nd part
	bool = xpl <= t[6] & xpl > t[5]
	b = (t[j+3] - xpl[bool])**2 / (t[j+3] - t[j+1]) / (t[j+3] - t[j+2])
	bpl[bool, 2] = (t[i+4] - xpl[bool]) / (t[i+4] - t[i+1]) * b

#3rd edge 1st part
	i = 3; j = i+1
	bool = xpl <= t[5] & xpl > t[4]
	a = (xpl[bool] - t[i]) * (t[i+2] - xpl[bool]) / (t[i+2] - t[i]) / (t[i+2] - t[i+1]) + (t[i+3] - xpl[bool]) * (xpl[bool] - t[i+1]) / (t[i+3] - t[i+1]) / (t[i+2] - t[i+1])
	b = (xpl[bool] - t[j])**2 / (t[j+2] - t[j]) / (t[j+1] - t[j])
	bpl[bool, 3] = (xpl[bool] - t[i]) / (t[i+3] - t[i]) * a + (t[i+4] - xpl[bool]) / (t[i+4] - t[i+1]) * b

#3rd edge 2nd part
	bool = xpl <= t[6] & xpl > t[5]
	a = (t[i+3] - xpl[bool])**2 / (t[i+3] - t[i+1]) / (t[i+3] - t[i+2])
	b =  (xpl[bool] - t[j]) * (t[j+2] - xpl[bool]) / (t[j+2] - t[j]) / (t[j+2] - t[j+1]) + (t[j+3] - xpl[bool]) * (xpl[bool] - t[j+1]) / (t[j+3] - t[j+1]) / (t[j+2] - t[j+1])
	bpl[bool, 3] = (xpl[bool] - t[i]) / (t[i+3] - t[i]) * a + (t[i+4] - xpl[bool]) / (t[i+4] - t[i+1]) * b

#3rd edge 3rd part
	bool = xpl <= t[7] & xpl > t[6]
	b = (t[j+3] - xpl[bool])**2 / (t[j+3] - t[j+1]) / (t[j+3] - t[j+2])
	bpl[bool, 3] = (t[i+4] - xpl[bool]) / (t[i+4] - t[i+1]) * b

#4th to (m-1)th edge
	for (i in 4:(m-1)) {
		j = i+1
#1st part	
		bool = xpl <= t[i+1] & xpl > t[i]
		a =  (xpl[bool] - t[i])**2 / (t[i+2] - t[i]) / (t[i+1] - t[i])
		bpl[bool, i] = (xpl[bool] - t[i]) / (t[i+3] - t[i]) * a

#2nd part
		bool = xpl <= t[i+2] & xpl > t[i+1]
		a = (xpl[bool] - t[i]) * (t[i+2] - xpl[bool]) / (t[i+2] - t[i]) / (t[i+2] - t[i+1]) + (t[i+3] - xpl[bool]) * (xpl[bool] - t[i+1]) / (t[i+3] - t[i+1]) / (t[i+2] - t[i+1])
		b = (xpl[bool] - t[j])**2 / (t[j+2] - t[j]) / (t[j+1] - t[j])
		bpl[bool, i] = (xpl[bool] - t[i]) / (t[i+3] - t[i]) * a + (t[i+4] - xpl[bool]) / (t[i+4] - t[i+1]) * b

#3rd part
		bool = xpl <= t[i+3] & xpl > t[i+2]
		a =  (t[i+3] - xpl[bool])**2 / (t[i+3] - t[i+1]) / (t[i+3] - t[i+2])
		b = (xpl[bool] - t[j]) * (t[j+2] - xpl[bool]) / (t[j+2] - t[j]) / (t[j+2] - t[j+1]) + (t[j+3] - xpl[bool]) * (xpl[bool] - t[j+1]) / (t[j+3] - t[j+1]) / (t[j+2] - t[j+1])
		bpl[bool, i] = (xpl[bool] - t[i]) / (t[i+3] - t[i]) * a + (t[i+4] - xpl[bool]) / (t[i+4] - t[i+1]) * b

#4th part
		bool = xpl <= t[i+4] & xpl > t[i+3]
		b =  (t[j+3] - xpl[bool])**2 / (t[j+3] - t[j+1]) / (t[j+3] - t[j+2])
		bpl[bool, i] = (t[i+4] - xpl[bool]) / (t[i+4] - t[i+1]) * b
}

#mth edge 1st part
	i = m; j = i+1
 	bool = xpl <= t[m+1] & xpl > t[m]
	a = (xpl[bool] - t[i])**2 / (t[i+2] - t[i]) / (t[i+1] - t[i])
	bpl[bool, i] = (xpl[bool] - t[i]) / (t[i+3] - t[i]) * a

#mth edge 2nd part
	bool = xpl <= t[m+2] & xpl > t[m+1]
	a = (xpl[bool] - t[i]) * (t[i+2] - xpl[bool]) / (t[i+2] - t[i]) / (t[i+2] - t[i+1]) + (t[i+3] - xpl[bool]) * (xpl[bool] - t[i+1]) / (t[i+3] - t[i+1]) / (t[i+2] - t[i+1])
	b = (xpl[bool] - t[j])**2 / (t[j+2] - t[j]) / (t[j+1] - t[j])
	bpl[bool, i] = (xpl[bool] - t[i]) / (t[i+3] - t[i]) * a + (t[i+4] - xpl[bool]) / (t[i+4] - t[i+1]) * b

#mth edge 3rd part
	bool = xpl <= t[m+3] & xpl > t[m+2]
	a = (t[i+3] - xpl[bool])**2 / (t[i+3] - t[i+1]) / (t[i+3] - t[i+2])
	b = (xpl[bool] - t[j]) * (t[j+2] - xpl[bool]) / (t[j+2] - t[j]) / (t[j+2] - t[j+1]) + (t[j+3] - xpl[bool]) * (xpl[bool] - t[j+1]) / (t[j+3] - t[j+1]) / (t[j+2] - t[j+1])
	bpl[bool, i] = (xpl[bool] - t[i]) / (t[i+3] - t[i]) * a + (t[i+4] - xpl[bool]) / (t[i+4] - t[i+1]) * b

#(m+1)th edge 1st part
	i = m+1; j = i+1
	bool = xpl <= t[m+2] & xpl > t[m+1]
	a = (xpl[bool] - t[i])**2 / (t[i+2] - t[i]) / (t[i+1] - t[i])
	bpl[bool, i] = (xpl[bool] - t[i]) / (t[i+3] - t[i]) * a 

#(m+1)th edge 2nd part
	bool = xpl <= t[m+3] & xpl > t[m+2]
	a = (xpl[bool] - t[i]) * (t[i+2] - xpl[bool]) / (t[i+2] - t[i]) / (t[i+2] - t[i+1]) + (t[i+3] - xpl[bool]) * (xpl[bool] - t[i+1]) / (t[i+3] - t[i+1]) / (t[i+2] - t[i+1])
	b = (xpl[bool] - t[j])**2 / (t[j+2] - t[j]) / (t[j+1] - t[j])
	bpl[bool, i] = (xpl[bool] - t[i]) / (t[i+3] - t[i]) * a + (t[i+4] - xpl[bool]) / (t[i+4] - t[i+1]) * b

#(m+2)th edge
	i = m+2
	bool = xpl <= t[m+3] & xpl > t[m+2]
	a = (xpl[bool] - t[i])**2 / (t[i+2] - t[i]) / (t[i+1] - t[i])
	bpl[bool, m+2] = (xpl[bool] - t[i]) / (t[i+3] - t[i]) * a
}
	ans = new.env()
	ans$bmat = bmat
	ans$bpl = bpl
	ans$xpl = xpl
	ans$knots = knots
	ans$secder = secder 
	ans

}

###########################
#bcspline first derivative#
###########################
bcsplfirderi = function(x, knots = NULL, xmin = 0, xmax = 1) {
	m = length(knots) 
	t = 1:(m + 6)*0
	t[1:3] = xmin #min(x)
	t[(m+4):(m+6)] = xmax #max(x)
	t[4:(m+3)] = knots

  	nx = length(x)
	firder =  matrix(0, nrow = nx, ncol = (m + 2)) 

#1st edge: t4 and t5 
	i = 1; j = i+1
	bool = x < t[5] & x >= t[4]
	firder[bool, 1] = -3 * (t[5] - x[bool])^2 / (t[5] - t[2]) / (t[j+3] - t[j+1]) / (t[j+3] - t[j+2])

#2nd edge: t4, t5 and t6 
	i = 2; j = i+1
	bool = x <= t[5] & x >= t[4]
	#if (bool) { 
     	a = (t[i+3] - x[bool])**2 / (t[i+3] - t[i+1]) / (t[i+3] - t[i+2])
	b = (x[bool] - t[j]) * (t[j+2] - x[bool]) / (t[j+2] - t[j]) / (t[j+2] - t[j+1]) + (t[j+3] - x[bool]) * (x[bool] - t[j+1]) / (t[j+3] - t[j+1]) / (t[j+2] - t[j+1])
     	a1 = -2 * (t[i+3] - x[bool]) / (t[i+3] - t[i+1]) / (t[i+3] - t[i+2])
     	b1 = ((t[j+2] - x[bool])- (x[bool] - t[j])) / (t[j+2] - t[j]) / (t[j+2] - t[j+1]) + ((t[j+3] - x[bool]) - (x[bool] - t[j+1])) / (t[j+3] - t[j+1]) / (t[j+2] - t[j+1]) 
	firder[bool, 2] = a / (t[i+3] - t[i]) + a1 * (x[bool] - t[i]) / (t[i+3] - t[i]) - b / (t[i+4] - t[i+1]) + (t[i+4] - x[bool]) / (t[i+4] - t[i+1]) * b1
	#}	
	bool = x <= t[6] & x > t[5]
	#if (bool) {
	b = (t[j+3] - x[bool])**2 / (t[j+3] - t[j+1]) / (t[j+3] - t[j+2])
     	b1 =  -2 * (t[j+3] - x[bool]) / (t[j+3] - t[j+1]) / (t[j+3] - t[j+2])
     	firder[bool, 2] = -b / (t[i+4] - t[i+1]) + (t[i+4] - x[bool]) / (t[i+4] - t[i+1]) * b1
	#}

#3rd edge: t4, t5, t6 and t7 
	i = 3; j = i+1
	bool =  x <= t[5] & x >= t[4]
	#if (bool) {
	a = (x[bool] - t[i]) * (t[i+2] - x[bool]) / (t[i+2] - t[i]) / (t[i+2] - t[i+1]) + (t[i+3] - x[bool]) * (x[bool] - t[i+1]) / (t[i+3] - t[i+1]) / (t[i+2] - t[i+1])
	b = (x[bool] - t[j])**2 / (t[j+2] - t[j]) / (t[j+1] - t[j]) 
     	a1 = ((t[i+2] - x[bool]) - (x[bool] - t[i])) / (t[i+2] - t[i]) / (t[i+2] - t[i+1]) + ((t[i+3] - x[bool]) - (x[bool] - t[i+1])) / (t[i+3] - t[i+1]) / (t[i+2] - t[i+1])
     	b1 = 2 * (x[bool] - t[j]) / (t[j+2] - t[j]) / (t[j+1] - t[j])  
	firder[bool, 3] = a / (t[i+3] - t[i]) + (x[bool] - t[i]) / (t[i+3] - t[i]) * a1 - b / (t[i+4] - t[i+1]) + (t[i+4] - x[bool]) / (t[i+4] - t[i+1]) * b1 
	#}	
	bool = x <= t[6] & x > t[5]
	#if (bool) {
	a = (t[i+3] - x[bool])**2 / (t[i+3] - t[i+1]) / (t[i+3] - t[i+2])
	b = (x[bool] - t[j]) * (t[j+2] - x[bool]) / (t[j+2] - t[j]) / (t[j+2] - t[j+1]) + (t[j+3] - x[bool]) * (x[bool] - t[j+1]) / (t[j+3] - t[j+1]) / (t[j+2] - t[j+1])
     	a1 =  -2 * (t[i+3] - x[bool]) / (t[i+3] - t[i+1]) / (t[i+3] - t[i+2])
     	b1 = ((t[j+2] - x[bool]) - (x[bool] - t[j])) / (t[j+2] - t[j]) / (t[j+2] - t[j+1]) + ((t[j+3] - x[bool]) - (x[bool] - t[j+1])) / (t[j+3] - t[j+1]) / (t[j+2] - t[j+1])   
     	firder[bool, 3] =  a / (t[i+3] - t[i]) + (x[bool] - t[i]) / (t[i+3] - t[i]) * a1 - b / (t[i+4] - t[i+1]) + (t[i+4] - x[bool]) / (t[i+4] - t[i+1]) * b1 
	#}
	bool = x <= t[7] & x > t[6]
	#if (bool) {
    	b = (t[j+3] - x[bool])**2 / (t[j+3] - t[j+1]) / (t[j+3] - t[j+2])
    	b1 = -2 *(t[j+3] - x[bool]) / (t[j+3] - t[j+1]) / (t[j+3] - t[j+2])
	firder[bool, 3] = -b / (t[i+4] - t[i+1]) + (t[i+4] - x[bool]) / (t[i+4] - t[i+1]) * b1
	#}

#4th  to (m-1)th edge: t[i], t[i+1], t[i+2], t[i+3] and t[i+4] 
if (m > 4) {
	for (i in 4:(m-1)) {
    		j = i+1
		bool = x <= t[i+1] & x >= t[i]
		#if (bool) {
		a = (x[bool] - t[i])**2 / (t[i+2] - t[i]) / (t[i+1] - t[i])
       		a1 = 2 * (x[bool] - t[i]) / (t[i+2] - t[i]) / (t[i+1] - t[i])
		firder[bool, i] = a / (t[i+3] - t[i]) + (x[bool] - t[i]) / (t[i+3] - t[i]) * a1
		#}	
		bool = x <= t[i+2] & x > t[i+1]
		#if (bool) {
		a = (x[bool] - t[i]) * (t[i+2] - x[bool]) / (t[i+2] - t[i]) / (t[i+2] - t[i+1]) + (t[i+3] - x[bool]) * (x[bool] - t[i+1]) / (t[i+3] - t[i+1]) / (t[i+2] - t[i+1])
		b = (x[bool] - t[j])**2 / (t[j+2] - t[j]) / (t[j+1] - t[j])
       		a1 =  ((t[i+2] - x[bool]) - (x[bool] - t[i])) / (t[i+2] - t[i]) / (t[i+2] - t[i+1]) + ((t[i+3] - x[bool]) - (x[bool] - t[i+1])) / (t[i+3] - t[i+1]) / (t[i+2] - t[i+1])
       		b1 = 2 * (x[bool] - t[j]) / (t[j+2] - t[j]) / (t[j+1] - t[j])
		firder[bool, i] =  a / (t[i+3] - t[i]) + (x[bool] - t[i]) / (t[i+3] - t[i]) * a1 - b / (t[i+4] - t[i+1]) + (t[i+4] - x[bool]) / (t[i+4] - t[i+1]) * b1
		#}	
		bool = x <= t[i+3] & x > t[i+2] 
		#if (bool) {
     		a =  (t[i+3] - x[bool])**2 / (t[i+3] - t[i+1]) / (t[i+3] - t[i+2])
 	    	b = (x[bool] - t[j]) * (t[j+2] - x[bool]) / (t[j+2] - t[j]) / (t[j+2] - t[j+1]) + (t[j+3] - x[bool]) * (x[bool] - t[j+1]) / (t[j+3] - t[j+1]) / (t[j+2] - t[j+1])
      		a1 = -2 * (t[i+3] - x[bool]) / (t[i+3] - t[i+1]) / (t[i+3] - t[i+2])
      		b1 = ((t[j+2] - x[bool]) - (x[bool] - t[j])) / (t[j+2] - t[j]) / (t[j+2] - t[j+1]) + ((t[j+3] - x[bool]) - (x[bool] - t[j+1])) / (t[j+3] - t[j+1]) / (t[j+2] - t[j+1])
		firder[bool, i] = a / (t[i+3] - t[i]) + (x[bool] - t[i]) / (t[i+3] - t[i]) * a1 - b / (t[i+4] - t[i+1]) + (t[i+4] - x[bool]) / (t[i+4] - t[i+1]) * b1
		#}
		bool = x <= t[i+4] & x > t[i+3]
		#if (bool) {
    		b = (t[j+3] - x[bool])**2 / (t[j+3] - t[j+1]) / (t[j+3] - t[j+2])
      		b1 = -2 * (t[j+3] - x[bool]) / (t[j+3] - t[j+1]) / (t[j+3] - t[j+2])
		firder[bool, i] =  -b / (t[i+4] - t[i+1]) + (t[i+4] - x[bool]) / (t[i+4] - t[i+1]) * b1
		#}
	}
}
#mth edge: t[m], t[m+1], t[m+2] and t[m+3] #checked
  	i = m; j = i+1
	bool = x <= t[m+1] & x >= t[m]
	#if (bool) {
    	a = (x[bool] - t[i])**2 / (t[i+2] - t[i]) / (t[i+1] - t[i])
    	a1 = 2 * (x[bool] - t[i]) / (t[i+2] - t[i]) / (t[i+1] - t[i])
	firder[bool, m] = a / (t[i+3] - t[i]) + (x[bool] - t[i]) / (t[i+3] - t[i]) * a1
	#}	
	bool = x <= t[m+2] & x > t[m+1]
	#if (bool) {
	a = (x[bool] - t[i]) * (t[i+2] - x[bool]) / (t[i+2] - t[i]) / (t[i+2] - t[i+1]) + (t[i+3] - x[bool]) * (x[bool] - t[i+1]) / (t[i+3] - t[i+1]) / (t[i+2] - t[i+1])
	b = (x[bool] - t[j])**2 / (t[j+2] - t[j]) / (t[j+1] - t[j])
    	a1 =  ((t[i+2] - x[bool]) - (x[bool] - t[i])) / (t[i+2] - t[i]) / (t[i+2] - t[i+1]) + ((t[i+3] - x[bool]) - (x[bool] - t[i+1]) ) / (t[i+3] - t[i+1]) / (t[i+2] - t[i+1]) 
    	b1 = 2 * (x[bool] - t[j]) / (t[j+2] - t[j]) / (t[j+1] - t[j])
	firder[bool, m] = a / (t[i+3] - t[i]) + (x[bool] - t[i]) / (t[i+3] - t[i]) * a1 - b / (t[i+4] - t[i+1]) + (t[i+4] - x[bool]) / (t[i+4] - t[i+1]) * b1
	#}	
	bool = x <= t[m+3] & x > t[m+2] 
	#if (bool) {
    	a = (t[i+3] - x[bool])**2 / (t[i+3] - t[i+1]) / (t[i+3] - t[i+2])
	b = (x[bool] - t[j]) * (t[j+2] - x[bool]) / (t[j+2] - t[j]) / (t[j+2] - t[j+1]) + (t[j+3] - x[bool]) * (x[bool] - t[j+1]) / (t[j+3] - t[j+1]) / (t[j+2] - t[j+1])
    	a1 =  -2 * (t[i+3] - x[bool]) / (t[i+3] - t[i+1]) / (t[i+3] - t[i+2])
    	b1 = ((t[j+2] - x[bool]) - (x[bool] - t[j])) / (t[j+2] - t[j]) / (t[j+2] - t[j+1]) + ((t[j+3] - x[bool]) - (x[bool] - t[j+1])) / (t[j+3] - t[j+1]) / (t[j+2] - t[j+1])
	firder[bool, m] = a / (t[i+3] - t[i]) + (x[bool] - t[i]) / (t[i+3] - t[i]) * a1  - b / (t[i+4] - t[i+1]) + (t[i+4] - x[bool]) / (t[i+4] - t[i+1]) * b1  
	#}	

#(m+1)th edge: t[m+1], t[m+2] and t[m+3] 
  	i = m+1; j = i+1
	bool = x <= t[m+2] & x >= t[m+1]
	#if (bool) {
    	a = (x[bool] - t[i])**2 / (t[i+2] - t[i]) / (t[i+1] - t[i])
    	a1 = 2 * (x[bool] - t[i]) / (t[i+2] - t[i]) / (t[i+1] - t[i])
	firder[bool, m+1] = a / (t[i+3] - t[i]) + (x[bool] - t[i]) / (t[i+3] - t[i]) * a1
	#}
	bool = x <= t[m+3] & x > t[m+2]
	#if (bool) {
    	a = (x[bool] - t[i]) * (t[i+2] - x[bool]) / (t[i+2] - t[i]) / (t[i+2] - t[i+1]) + (t[i+3] - x[bool]) * (x[bool] - t[i+1]) / (t[i+3] - t[i+1]) / (t[i+2] - t[i+1])
	b = (x[bool] - t[j])**2 / (t[j+2] - t[j]) / (t[j+1] - t[j])
    	a1 = ((t[i+2] - x[bool]) - (x[bool] - t[i])) / (t[i+2] - t[i]) / (t[i+2] - t[i+1]) + ((t[i+3] - x[bool]) - (x[bool] - t[i+1])) / (t[i+3] - t[i+1]) / (t[i+2] - t[i+1])
    	b1 = 2 * (x[bool] - t[j]) / (t[j+2] - t[j]) / (t[j+1] - t[j])
	firder[bool, m+1] = a / (t[i+3] - t[i]) + (x[bool] - t[i]) / (t[i+3] - t[i]) * a1 - b / (t[i+4] - t[i+1]) + (t[i+4] - x[bool]) / (t[i+4] - t[i+1]) * b1  
	#}

#(m+2)th edge: t[m+2] and t[m+3]
  	i = m+2
	bool = x <= t[m+3] & x >= t[m+2]
	#if (bool) {
    	a = (x[bool] - t[i])**2 / (t[i+2] - t[i]) / (t[i+1] - t[i])
    	a1 = 2 * (x[bool] - t[i]) / (t[i+2] - t[i]) / (t[i+1] - t[i]) 
	firder[bool, m+2] =  a / (t[i+3] - t[i]) + (x[bool] - t[i]) / (t[i+3] - t[i]) * a1
	#}	
	firder
}

###########
#CicFamily#
###########
CicFamily <- function(object,...)UseMethod("CicFamily")
CicFamily <- function(object) {
  llh.fun <- function(y, muhat = NULL, etahat = NULL, n = NULL, weights = NULL, fml = object$family){
    sm <- 1e-7
    #sm <- 1e-5
    if (is.null(weights)) {
	weights <- 1:n*0 + 1
    }
    w <- weights
    if (fml == "poisson") {
      llh <- 2 * sum(w * (muhat - y * etahat)) / n
    }
    if (fml == "binomial") {
      llh <- 0
      if (all(0 <= y) & all(y <= 1)) {
        for (i in 1:n) {
          if (muhat[i] > 0 & muhat[i] < 1) {
            llh <- llh + w[i] * (y[i] * log(muhat[i]) + (1 - y[i]) * log(1 - muhat[i])) 
          }
        }
        llh <- (-2/n) * llh
      } else {
          stop ("y values must be 0 <= y <= 1!")
      }
    }
    if (fml == "gaussian") {
      if (all(w == 1)) {
        llh <- log(sum((y - etahat)^2))
      } else {
          llh <- log(sum(w * (y - etahat)^2)) - sum(log(w)) / n
      }
    }
    llh 
  }

  etahat.fun <- function(n, y, fml = object$family){
    if (fml == "poisson") {
      etahat <- 1:n*0 + log(mean(y)) 
    } 
    if (fml == "binomial") {
      etahat <- 1:n*0 
    }
    etahat
  }

  gr.fun <- function(y, etahat = NULL, weights = NULL, fml = object$family){
    n <- length(y)
    if (is.null(weights)) {
      weights <- 1:n*0 + 1
    }
    w <- weights 
    if (fml == "poisson") {
       gr <- w * (exp(etahat) -  y) 
    }
    if (fml == "binomial") {
       if (all(etahat == 0)) { 
         gr <- w * (1/2 - y)
       } else {
	   gr <- 1:n*0
	   for (i in 1:n) {
	     if (etahat[i] > 100) {
		 gr[i] <- w[i] * (1 - y[i]) 
	     } else {gr[i] <- w[i] * (exp(etahat[i]) / (1 + exp(etahat[i])) - y[i])}
           }         
         }
    }
    gr
  }

  wt.fun <- function(etahat = NULL, n = NULL, weights = NULL, fml = object$family){
    if (is.null(weights)) {
	weights <- 1:n*0 + 1
    }
    w <- weights 
    if (fml == "poisson") {
      wt <-  w * exp(etahat)
    }
    if (fml == "binomial") {
      if (all(etahat == 0)){
        #wt <- 1:n*0 + 1/4
	wt <- w * (1:n*0 + 1/4)
      } else {
	  wt <- 1:n*0
          for (i in 1:n) {
            if (etahat[i] > 100) {
              wt[i] <- 0
            } else {
                wt[i] <- w[i] * exp(etahat[i]) / ((1 + exp(etahat[i]))^2)
              }
          }
        }
    }
    if (fml == "gaussian") {
      wt <- w # (1:n*0 + 1) / w 
    }
    wt <- as.vector(wt)
    wt 
  }

  zvec.fun <- function(cvec = NULL, wt = NULL, y, sm = 1e-7, fml = object$family) {
    n <- length(y)
    if (fml == "gaussian") {
      #zvec <- y
      zvec <- wt^(1/2) * y
    }
    if (fml == "poisson") {	
      #zvec <- cvec / wt
      zvec <- cvec / sqrt(wt) 
    }
    if (fml == "binomial") {
     zvec = 1:n*0
     zvec[wt == 0] <- 1 / sm
     zvec[wt > 0] <- cvec[wt > 0] / sqrt(wt[wt > 0])
    }
    zvec 
  }

  muhat.fun <- function(etahat, wt = NULL, fml = object$family){
    n <- length(etahat)
    if (fml == "poisson") {
      muhat <- exp(etahat)
    }
    if (fml == "binomial") {
      muhat <- 1:n*0
      #muhat[wt == 0] <- 1
      for (i in 1:n) {
        if (etahat[i] > 100) {
	  muhat[i] <- 1
 	} else {
          muhat[i] <- exp(etahat[i]) / (1 + exp(etahat[i]))
        }
      }
    }
    if (fml == "gaussian") {
      muhat <- etahat
    }
   muhat 
  }

  ysim.fun <- function(n, mu0 = NULL, fml = object$family) {
    if (fml == "binomial") {
      ysim <- 1:n*0
      ysim[runif(n) < .5] <- 1
    }
    if (fml == "poisson") {
      if (!is.null(mu0)) {
        ysim <- rpois(n, mu0)
      }
    }
    if (fml == "gaussian") {
      ysim <- rnorm(n)
    }
    ysim 
  }

  deriv.fun <- function(muhat, fml = object$family) {
    if (fml == "binomial") {
	deriv <- 1 / (muhat * (1 - muhat))
    }
    if (fml == "poisson") {
	deriv <- 1 / muhat
    }
    if (fml == "gaussian") {
	deriv <- 1
    }
   deriv
  }

 dev.fun <- function(y, muhat, etahat, weights, fml = object$family){
  n <- length(y)
  sm <- 1e-7
  #sm <- 1e-5
  if (is.null(weights)) {
	weights <- 1:n*0 + 1
  }
  w <- weights
  vmat <- matrix(1:n*0 + 1, ncol = 1)
  if (fml == "poisson") {
        #dev <- 2 * sum(w * (y * log(y / muhat) - y + muhat))
	dev <- 0
	for (i in 1:n) {
	  if (y[i] == 0) {
            dev <- dev + 2 * w[i] * muhat[i]
          } else {
            dev <- dev + 2 * w[i] * (y[i] * log(y[i] / muhat[i]) - y[i] + muhat[i])
          }
	}
  }
  if (fml == "binomial") {
        dev <- 0
        for (i in 1:n) {
          if (y[i] == 0) {
            dev <- dev + 2 * w[i] * log(w[i] / (w[i] - w[i] * muhat[i]))
          } else if (y[i] == 1) {
              dev <- dev + 2 * w[i] * log(w[i] / (w[i] * muhat[i]))
          } else if (0 < y[i] & y[i] < 1) {
              dev <- dev + 2 * w[i] * y[i] * log(w[i] * y[i] / (w[i] * muhat[i])) + 2 * (w[i] - w[i] * y[i]) * log((w[i] - w[i] * y[i]) / (w[i] - w[i] * muhat[i]))
          } else {
             stop ("y values must be 0 <= y <= 1!")
          }
       }
  }
  if (fml == "gaussian") {
        dev <- sum(w * (y - muhat)^2)
  }
###################
#get null deviance#
###################
  if (fml == "binomial" | fml == "poisson") {
      diff <- 1
      muhat0 <- mean(y) + 1:n*0
      if (fml == "poisson") {
         etahat0 <- log(muhat0)
      } 
      if (fml == "binomial") {
         etahat0 <- log(muhat0 / (1 - muhat0))
      } 		
      while (diff > sm) {
        oldmu <- muhat0
	zhat <- etahat0 + (y - muhat0) * deriv.fun(muhat0, fml = fml)		
	wmat <- diag(as.vector(w / deriv.fun(muhat0, fml = fml)))			
	b <- solve(t(vmat) %*% wmat %*% vmat) %*% t(vmat) %*% wmat %*% zhat
	etahat0 <- vmat %*% b
	muhat0 <- muhat.fun(etahat0, fml = fml)		
	diff <- mean((muhat0 - oldmu)^2)	
      }
      if (fml == "poisson") {
        #dev.null <- 2 * sum(w * (y * log(y / muhat0) - y + muhat0))
	dev.null <- 0
        for (i in 1:n) {
	  if (y[i] == 0) {
            dev.null <- dev.null + 2 * w[i] * muhat0[i]
          } else {
            dev.null <- dev.null + 2 * w[i] * (y[i] * log(y[i] / muhat0[i]) - y[i] + muhat0[i])
          }
	}
      }
      if (fml == "binomial") {
        dev.null <- 0
        for (i in 1:n) {
          if (y[i] == 0) {
            dev.null <- dev.null + 2 * w[i] * log(w[i] / (w[i] - w[i] * muhat0[i]))
          } else if (y[i] == 1) {
              dev.null <- dev.null + 2 * w[i] * log(w[i] / (w[i] * muhat0[i]))
          } else if (0 < y[i] & y[i] < 1) {
              dev.null <- dev.null + 2 * w[i] * y[i] * log(w[i] * y[i] / (w[i] * muhat0[i])) + 2 * (w[i] - w[i] * y[i]) * log((w[i] - w[i] * y[i]) / (w[i] - w[i] * muhat0[i]))
          } else {
              stop ("y values must be 0 <= y <= 1!")
	  }
        }
      } 
  }
  if (fml == "gaussian") {
     wmat <- diag(w)
     b <- solve(t(vmat) %*% wmat %*% vmat) %*% t(vmat) %*% wmat %*% y
     etahat0 <- vmat %*% b
     muhat0 <- muhat.fun(etahat0, fml = fml)	
     dev.null <- sum(w * (y - muhat0)^2)
  }
  rslt <- new.env()
  rslt$dev <- dev
  rslt$dev.null <- dev.null 
  rslt
  }
  ans <- list(llh.fun = llh.fun, etahat.fun = etahat.fun, gr.fun = gr.fun, wt.fun = wt.fun, zvec.fun = zvec.fun, muhat.fun = muhat.fun, ysim.fun = ysim.fun, deriv.fun = deriv.fun, dev.fun = dev.fun)
  class(ans) <- "CicFamily"
  return (ans)
}


##########
#make_arp#
##########
make_arp = function(y, pnt = FALSE, dmat = NULL, bmat = NULL, sl = NULL, p_max = 2, m2 = NULL, lambdas = NULL, pen_fit = NULL, p_fit = NULL, dv0 = NULL, qv0 = NULL, cv = NULL, hs = NULL, ids = NULL, constr = TRUE) {
	if (is.null(pen_fit) & is.null(p_fit)) {
		ps = 1:p_max
		if (!pnt) {
	   		lambdas = as.vector(0)
 		}
  		la = length(lambdas)
#include p = 0
  		nposs = la * (p_max + 1)
	} else {
  		ps = as.vector(p_fit)
  		lambdas = as.vector(pen_fit)
  		la = p_max = nposs = 1
	}
  	aics = aics1 = aics2 = trs = trrs = sigs = pens = 1:nposs*0
	n = length(y)
  	phis = bhs = r_mats = pmat0s = list()
	fhats = fhats_x = ess = rss = matrix(0, nrow = nposs, ncol = n)
	sm = 1e-5
	iter = 0
# for p: 0 ~ 4
  	rs = rs2 = 1:n*0
  	r_mat = tau_mat = matrix(0, n, n)
	qv0 = crossprod(bmat)
	tdf = ncol(bmat)
    	imat = diag(tdf)
  	for (ila in 1:la) {
        	pen = lambdas[ila]
         	if (pnt) {
            		qv = qv0 + pen * dv0
         	} else {
            		qv = qv0
         	}	
		cv = crossprod(bmat, y)
#new, for jp:
		if (constr) {
         		qans = qprog(qv, cv, sl, 1:nrow(sl)*0, msg = FALSE)	
			bh = fitted(qans)
         		qdf = qans$df
          		fhat0 = bmat %*% bh
          		d_use = t(qans$xmat)
		} else {
			bh = solve(qv0, cv)
			qdf = tdf 
			fhat0 = bmat %*% bh	
			d_use = bmat
		}          	
		if (qdf < tdf) {
             		pd = d_use %*% solve(crossprod(d_use), t(d_use))
             		pmat0 = imat - pd
          	} else {
            		pmat0 = imat
          	}                    
          	umat = chol(qv)
          	iumat = diag(ncol(umat))
          	uinv = backsolve(umat, iumat)
          	bu = bmat %*% uinv    
          	pmat = bu %*% tcrossprod(pmat0, bu)
#constrained
          	tr = sum(diag(pmat))
		pu = bmat %*% solve(qv, t(bmat))
#unconstrained
		trr = sum(diag(pu))  
          	es = y - fhat0
          	nd = ncol(d_use)
          	sig2 =  sum((y - fhat0)^2) / (n - 1.5 * (tdf - nd))
          	aic = n * log(sig2) + 2 * tr
  		if (is.null(p_fit) || p_fit == 0) {
   			iter = iter + 1
        		phis[[iter]] = 0
         		bhs[[iter]] = bh
			r_mats[[iter]] = diag(n)    		
			fhats[iter, ] = fhat0
         		ess[iter, ] = es
         		aics[iter] = aic		
         		trs[iter] = tr
			trrs[iter] = trr 
         		sigs[iter] = sig2
         		pens[iter] = pen    
		}
		if (is.null(p_fit) || p_fit != 0) {
			for (ip in 1:p_max) {
				nrep = 0
  				sm0 = 1
                		fhat = fhat0 
				iter = iter + 1  				
				p = ps[ip]  
  				while (sm0 > sm & nrep < 1e+3) {
                      			nrep = nrep + 1
                      			oldf = fhat 
                      			es = y - fhat
					es = es - mean(es)                      			
					rs = sapply(hs, FUN = function(h, X = es, N = n) {sum((X[1:(N - h)] - mean(X)) * (X[(1 + h):N]- mean(X))) / N})
					tau_mat = matrix(0, n, n)
                      			for (ih in 1:n) {
                        		  	ids_i = ids[[ih]]
	                        		h = hs[ih]
	                        		tau_mat[ids_i] = rs[h + 1]
                      			}               
					tau_p_mat = tau_mat[1:p, 1:p]
            				rp = rs[2:(p + 1)]
					phi = solve(tau_p_mat, rp)  
					sig2 = rs[1] - crossprod(phi, rp)        
					rs2 = make_rs(phi, ord = p, nvec = n, sig2 = sig2)
					r_mat = matrix(0, n, n)
     					for (ih in 1:n) {
                        		  	ids_i = ids[[ih]]
	                        		h = hs[ih]
	                        		r_mat[ids_i] = rs2[h + 1]
                      			}
#this line makes the fit flat!
					r_mat = r_mat / rs2[1]
					if (pnt) {
     	                   			qv = crossprod(bmat, solve(r_mat, bmat)) + pen * dv0
    		        		} else {
         					qv = crossprod(bmat, solve(r_mat, bmat))                        	
      		        		}
					cv = crossprod(bmat, solve(r_mat, y))
#new, for jp:
					if (constr) {
                      				qans = qprog(qv, cv, sl, 1:nrow(sl)*0, msg = FALSE)
                      				bh = fitted(qans)
                      				fhat = bmat %*% bh
					} else {
						bh = solve(qv, cv)
						fhat = bmat %*% bh	
					}
#lines(x, fhat, lty = 2, col = nrep)
                      			sm0 = mean((fhat - oldf)^2)	
            			}
				if (constr) {		
					d_use = t(qans$xmat)
            				qdf = qans$df
				} else {
					d_use = bmat
					qdf = tdf 
				}
            			if (qdf < tdf) {
					pd = d_use %*% solve(crossprod(d_use), t(d_use))
               				pmat0 = imat - pd
            			} else {
              				pmat0 = imat
            			}   
				umat = chol(qv)
				iumat = diag(ncol(umat))
            			uinv = backsolve(umat, iumat)
           			bu = bmat %*% uinv
            			pmat = bu %*% tcrossprod(pmat0, bu) %*% chol2inv(chol(r_mat))    
				pu = bmat %*% solve(qv, crossprod(bmat, chol2inv(chol(r_mat))))
#constrained
            			tr = sum(diag(pmat))
#unconstrained	
				trr = sum(diag(pu))   			
				rss[iter, ] = rs2
				pmat0s[[iter]] = pmat0         		
         			bhs[[iter]] = bh
	 			fhats[iter, ] = fhat              
         			ess[iter, ] = es 
        			r_mats[[iter]] = r_mat
        			phis[[iter]] = phi
	 			aics1[iter] = n * log(sig2)
            			aics2[iter] = 2 * (p + tr) 	
            			aics[iter] = n * log(sig2) + 2 * (p + tr)	
         			trs[iter] = tr
				trrs[iter] = trr
         			sigs[iter] = sig2
         			pens[iter] = pen
   			}
     		}                     
	}
  	id_aic = which(aics == min(aics))
  	lid = length(id_aic)
  	if (lid > 1) {
     		id_use = id_aic[1]
  	} else {
     		id_use = id_aic
  	}
#print (id_use)
	phi_use = phis[[id_use]]
	r_use = r_mats[[id_use]]
  	bh_use = bhs[[id_use]]
	fhat_use = fhats[id_use, ]
  	tr_use = trs[id_use]
	trr_use = trrs[id_use]
  	sig_use = sigs[id_use]
  	es_use = ess[id_use, ]
  	lambda_use = pens[id_use]
#lambda by order
	sigs = matrix(sigs, nrow = la, byrow = TRUE)
	aics = matrix(aics, nrow = la, byrow = TRUE)
	trs = matrix(trs, nrow = la, byrow = TRUE)
	trrs = matrix(trrs, nrow = la, byrow = TRUE)
	pens = matrix(pens, nrow = la, byrow = TRUE) 	
	aicmin = min(aics)
  	rslt = list(phi_use = phi_use, phis = phis, thetahat = bh_use, fhat_use = fhat_use, fhats = fhats, df = tr_use, r_use = r_use, sig_use = sig_use, sigs = sigs, pens = pens, trs = trs, trrs = trrs, aics = aics, aics1 = aics1, aics2 = aics2, aicmin = aicmin, tr_use = tr_use, trr_use = trr_use, es_use = es_use, ess = ess, lambda_use = lambda_use)
  	rslt
}


############################
#make theoretical residuals#
############################
make_rs = function(phi, ord = 1, nvec, sig2) {
	rhs = rvec = 1:nvec*0
	if (ord == 1) {
		r0 = sig2 / (1 - phi^2)
		rhs[1] = 1
		rhs[2:nvec] = phi^(1:(nvec - 1))
		rvec = rhs * r0	
		#rvec = rhs	
	}
	if (ord == 2) {
		phi1 = phi[1]; phi2 = phi[2]
		rh1 = phi1 / (1 - phi2); rh2 = phi1^2 / (1 - phi2) + phi2
		rhs[1] = 1; rhs[2] = rh1; rhs[3] = rh2
		rh_vec = c(rh1, rh2)
		r0 = sig2 / (1 - sum(phi * rh_vec))
		if (nvec >= 3) {
			rh_iter = rev(rh_vec)
			for (l in 1:(nvec - 3)) {
				rhs[3 + l] = sum(phi * rh_iter)
				rh_iter = c(rhs[3 + l], rh_iter[1])
				if (all(abs(rh_iter) < 1e-10)) {			
					break
				}
			}	
		}
		rvec = rhs * r0
		#rvec = rhs	
	}
	if (ord == 3) {
		phi_mat = matrix(0, 3, 3)
		phi_mat[1, 1] = 1 - phi[2]; phi_mat[1, 2] = -phi[3]
		phi_mat[2, 1] = -phi[1] - phi[3]; phi_mat[2, 2] = 1
		phi_mat[3, 1] = -phi[2]; phi_mat[3, 2] = -phi[1]; phi_mat[3, 3] = 1
		rh_vec = solve(phi_mat, phi)
		rhs[1] = 1
		rhs[2:4] = rh_vec
		r0 = sig2 / (1 - sum(phi * rh_vec))
		if (nvec >= 4) {
			rh_iter = rev(rh_vec)
			for (l in 1:(nvec - 4)) {
				rhs[4 + l] = sum(phi * rh_iter)
				rh_iter = c(rhs[4 + l], rh_iter[1:2])
				if (all(abs(rh_iter) < 1e-10)) {
					break
				}
			}	
		}
		rvec = rhs * r0
		#rvec = rhs	
	}
	if (ord == 4) {
		phi_mat = matrix(0, 4, 4)
		phi_mat[1, 1] = 1 - phi[2]; phi_mat[1, 2] = -phi[3]; phi_mat[1, 3] = -phi[4]
		phi_mat[2, 1] = -phi[1] - phi[3]; phi_mat[2, 2] = -phi[4] + 1
		phi_mat[3, 1] = -phi[2] - phi[4]; phi_mat[3, 2] = -phi[1]; phi_mat[3, 3] = 1
		phi_mat[4, 1] = -phi[3]; phi_mat[4, 2] = -phi[2]; phi_mat[4, 3] = -phi[1]; phi_mat[4, 4] = 1
		rh_vec = solve(phi_mat, phi)
		rhs[1] = 1
		rhs[2:5] = rh_vec
		r0 = sig2 / (1 - sum(phi * rh_vec))
		if (nvec >= 5) {
			rh_iter = rev(rh_vec)
			for (l in 1:(nvec - 5)) {
				rhs[5 + l] = sum(phi * rh_iter)
				rh_iter = c(rhs[5 + l], rh_iter[1:3]) 
				if (all(abs(rh_iter) < 1e-10)) {
					break
				}
			}	
		}
		rvec = rhs * r0
		#rvec = rhs
	}
	rvec = round(rvec, 8)
	return (rvec)
}

