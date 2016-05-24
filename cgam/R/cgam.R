######
#cgam#
######
cgam <- function(formula, nsim = 1e+2, family = gaussian(), data = NULL, weights = NULL)
{
  cl <- match.call()
  if (is.character(family)) 
     family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) 
     family <- family()
  if (is.null(family$family)) 
     stop("'family' not recognized!")
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  ynm <- names(mf)[1]
  mt <- attr(mf, "terms")
  y <- model.response(mf, "any")
  if (family$family == "binomial") {
	if (class(y) == "factor") {
		y = ifelse(y == levels(y)[1], 0, 1)
	}
  }
  shapes1 <- NULL; shapes2 <- NULL
  xmat <- NULL; xmatnms <- NULL
  tr <- NULL; umb <- NULL
  tree.delta <- NULL; umbrella.delta <- NULL 
  tid1 <- NULL; tid2 <- NULL; tpos2 <- 0
  uid1 <- NULL; uid2 <- NULL; upos2 <- 0
  nums <- NULL; ks <- list(); sps <- NULL; xid <- 1
  zmat <- NULL; zid <- NULL; zid0 <- NULL; zid1 <- NULL; zid2 <- NULL; znms <- NULL; is_param <- NULL; is_fac <- NULL; vals <- NULL; st <- 1; ed <- 1
  for (i in 2:ncol(mf)) {
    if (is.numeric(attributes(mf[,i])$shape)) {
       shapes1 <- c(shapes1, attributes(mf[,i])$shape)
       xmat <- cbind(xmat, mf[,i])
       xmatnms <- c(xmatnms, attributes(mf[,i])$nm)
       nums <- c(nums, attributes(mf[,i])$numknots)
       sps <- c(sps, attributes(mf[,i])$space)
       ks[[xid]] <- attributes(mf[,i])$knots
       xid <- xid + 1
    }
    if (is.character(attributes(mf[,i])$shape)) {
       shapes2 <- c(shapes2, attributes(mf[,i])$shape)
       if (attributes(mf[,i])$shape == "tree") {
		tree.delta <- rbind(tree.delta, tree.fun(mf[,i]))
		tpos1 <- tpos2 + 1
		tpos2 <- tpos2 + nrow(tree.fun(mf[,i]))
		tid1 <- c(tid1, tpos1)
		tid2 <- c(tid2, tpos2)
		tr <- cbind(tr, mf[,i])
       }
       if (attributes(mf[,i])$shape == "umbrella") {
		umbrella.delta <- rbind(umbrella.delta, umbrella.fun(mf[,i]))
		upos1 <- upos2 + 1
		upos2 <- upos2 + nrow(umbrella.fun(mf[,i]))
		uid1 <- c(uid1, upos1)
		uid2 <- c(uid2, upos2)
		umb <- cbind(umb, mf[,i])	
	}
    }
    if (is.null(attributes(mf[,i])$shape)) {
	if (!is.null(names(mf)[i])) {
	  znms <- c(znms, names(mf)[i])
	}
        if (!is.matrix(mf[,i])) {
          zid <- c(zid, i)
	  is_param <- c(is_param, TRUE)
          if (is.factor(mf[,i])) {
	    is_fac <- c(is_fac, TRUE)
	    ch_char <- suppressWarnings(is.na(as.numeric(levels(mf[, i]))))
            if (any(ch_char)) {
	      vals <- c(vals, unique(levels(mf[, i]))[-1])
            } else {
	      vals <- c(vals, as.numeric(levels(mf[, i]))[-1])
	    }
            nlvs <- length(attributes(mf[,i])$levels)
	    ed <- st + nlvs - 2 
	    zid1 <- c(zid1, st)
	    zid2 <- c(zid2, ed)
	    st <- st + nlvs - 1
	    zmat0 <- model.matrix(~ mf[, i])[, -1, drop = FALSE]
	    zmat <- cbind(zmat, zmat0)
          } else {
	    is_fac <- c(is_fac, FALSE)
            zmat <- cbind(zmat, mf[, i])
	    ed <- st
            zid1 <- c(zid1, st)
	    zid2 <- c(zid2, ed) 
	    st <- st + 1
	    vals <- c(vals, "")
         }
       } else {
	  is_param <- c(is_param, FALSE)
          is_fac <- c(is_fac, FALSE)
	  zmat0 <- mf[, i]
	  mat_cols <- ncol(zmat0)
	  mat_rm <- NULL
	  #rm_num <- 0
	  for (irm in 1:mat_cols) {
       	  	if (all(round(diff(zmat0[, irm]), 8) == 0)) {
                	mat_rm <- c(mat_rm, irm)
          	}
   	  }
	  if (!is.null(mat_rm)) {
	  	zmat0 <- zmat0[, -mat_rm, drop = FALSE]
		#rm_num <- rm_num + length(mat_rm)
	  }
	  zmat <- cbind(zmat, zmat0)
	  vals <- c(vals, 1)
	  zid <- c(zid, i)
	  nlvs <- ncol(zmat0) + 1
	  ed <- st + nlvs - 2
	  zid1 <- c(zid1, st)
	  zid2 <- c(zid2, ed)
	  st <- st + nlvs - 1
      }
    }
  }
  dimnames(zmat)[[2]] <- NULL
  if (family$family == "binomial" | family$family == "poisson") {
     wt.iter = TRUE
  } else {wt.iter = FALSE}
  if (is.null(shapes1) & is.null(shapes2)) {
    nsim <- 0
  }
  #attr(xmat, "shape") <- shapes1
  xmat0 <- xmat; shapes0 <- shapes1; nums0 <- nums; ks0 <- ks; sps0 <- sps; xmatnms0 <- xmatnms
  if (any(shapes1 == 17)) {
    kshapes <- length(shapes1)
    obs <- 1:kshapes
    idx_s <- obs[which(shapes1 == 17)]; idx <- obs[which(shapes1 != 17)]
  
    xmat0[ ,1:length(idx_s)] <- xmat[ ,idx_s]
    shapes0[1:length(idx_s)] <- shapes1[idx_s]
    nums0[1:length(idx_s)] <- nums[idx_s]
    sps0[1:length(idx_s)] <- sps[idx_s]
    ks0[1:length(idx_s)] <- ks[idx_s]
    xmatnms0[1:length(idx_s)] <- xmatnms[idx_s]

    if (length(idx) > 0) {
      xmat0[ ,(1 + length(idx_s)):kshapes] <- xmat[ ,idx]
      shapes0[(1 + length(idx_s)):kshapes] <- shapes1[idx]
      nums0[(1 + length(idx_s)):kshapes] <- nums[idx]
      sps0[(1 + length(idx_s)):kshapes] <- sps[idx]
      ks0[(1 + length(idx_s)):kshapes] <- ks[idx]
      xmatnms0[(1 +length(idx_s)):kshapes] <- xmatnms[idx]
    }
    xmat <- xmat0; nums <- nums0; ks <- ks0; sps <- sps0; xmatnms <- xmatnms0
  }
  shapes <- c(shapes1, shapes2)
  ans <- cgam.fit(y = y, xmat = xmat0, zmat = zmat, shapes = shapes0, numknots = nums0, knots = ks0, space = sps0, nsim = nsim, family = family, wt.iter = wt.iter, umbrella.delta = umbrella.delta, tree.delta = tree.delta, weights = weights)
  if (!is.null(uid1) & !is.null(uid2)) {
    uid1 <- uid1 + ans$d0 + ans$capm
    uid2 <- uid2 + ans$d0 + ans$capm
  }
  if (!is.null(tid1) & !is.null(tid2)) {
    tid1 <- tid1 + ans$d0 + ans$capm + ans$capu
    tid2 <- tid2 + ans$d0 + ans$capm + ans$capu 
  }
  rslt <- list(etahat = ans$etahat, muhat = ans$muhat, vcoefs = ans$vcoefs, xcoefs = ans$xcoefs, zcoefs = ans$zcoefs, ucoefs = ans$ucoefs, tcoefs = ans$tcoefs, coefs = ans$coefs, cic = ans$cic, d0 = ans$d0, edf0 = ans$edf0, etacomps = ans$etacomps, xmat = xmat, zmat = zmat, tr = tr, umb = umb, tree.delta = tree.delta, umbrella.delta = umbrella.delta, bigmat = ans$bigmat, shapes = shapes, wt = ans$wt, wt.iter = ans$wt.iter, family = ans$family, SSE0 = ans$sse0, SSE1 = ans$sse1, pvals.beta = ans$pvals.beta, se.beta = ans$se.beta, null_df = ans$df.null, df = ans$df, resid_df_obs = ans$resid_df_obs, null_deviance = ans$dev.null, deviance = ans$dev, tms = mt, capm = ans$capm, capms = ans$capms, capk = ans$capk, capt = ans$capt, capu = ans$capu, xid1 = ans$xid1, xid2 = ans$xid2, tid1 = tid1, tid2 = tid2, uid1 = uid1, uid2 = uid2, zid = zid, vals = vals, zid1 = zid1, zid2 = zid2, nsim = nsim, xmatnms = xmatnms,  ynm = ynm, znms = znms, is_param = is_param, is_fac = is_fac, knots = ans$knots, numknots = ans$numknots, sps = sps, ms = ans$ms)
  rslt$call <- cl
  class(rslt) <- "cgam"
  return (rslt) 
}

###############
#amat function#
###############
amat.fun <- function(x)
{
	obs <- 1:length(x)
	amat <- NULL
	xu <- unique(x)
	nx <- sort(xu)
	hd <- head(nx, 1)
	tl <- nx
	while (length(tl) > 1) {
		hd <- head(tl, 1)
		tl <- tl[-1]
		paired <- 0
		for (i in 1:length(tl)) {
			a1 <- 1:length(x)*0
			if (hd * tl[i] > 0) {
				if (hd < 0) {
		        		if (tl[i] > hd) {
						a1[min(obs[which(x == hd)])] <- -1; a1[min(obs[which(x == tl[i])])] <- 1 
					} else {
						a1[min(obs[which(x == hd)])] <- 1; a1[min(obs[which(x == tl[i])])] <- -1
					}
				}
				if (hd > 0) {
					if (tl[i] < hd) {
						a1[min(obs[which(x == hd)])] <- -1; a1[min(obs[which(x == tl[i])])] <- 1 
					} else {
						a1[min(obs[which(x == hd)])] <- 1; a1[min(obs[which(x == tl[i])])] <- -1
					}
				}
			#if (!all(a1 == 0) & paired == 0 ) {amat <- rbind(amat, a1); paired <- 1}
			}
			if (hd * tl[i] == 0) {
				if (hd == 0) {
					#if (tl[i] > 0) {
						a1[min(obs[which(x == hd)])] <- 1; a1[min(obs[which(x == tl[i])])] <- -1
					#}
				} else {
					a1[min(obs[which(x == hd)])] <- -1; a1[min(obs[which(x == tl[i])])] <- 1
				}
			}
			if (!all(a1 == 0) & paired == 0 ) {amat <- rbind(amat, a1); paired <- 1}
		}
	}
	dimnames(amat) <- NULL
	amat
}	

###############
#bmat function#
###############
bmat.fun <- function(x)
{
	obs <- 1:length(x)
	bmat <- NULL
	hd <- head(x, 1)
	tl <- x
	j <- 0
	while (length(tl) > 1) {
		hd <- head(tl, 1)
		tl <- tl[-1]
		paired <- 0
		j <- j + 1
		for (i in 1:length(tl)) {
			b1 <- 1:length(x)*0
			if (hd == tl[i]) {b1[j] <- -1; b1[j + i] <- 1}
			if (!all(b1 == 0) & paired == 0 ) {bmat <- rbind(bmat, b1); paired <- 1}
		}
	}
	dimnames(bmat) <- NULL
	bmat
}

##########
#cgam.fit#
##########
cgam.fit <- function(y, xmat, zmat, shapes, numknots, knots, space, nsim, family = gaussian(), wt.iter = FALSE, umbrella.delta = NULL, tree.delta = NULL, weights = NULL) {
#print (weights)
        linkfun <- family$linkfun
	cicfamily <- CicFamily(family)
	llh.fun <- cicfamily$llh.fun
	etahat.fun <- cicfamily$etahat.fun
	gr.fun <- cicfamily$gr.fun
	wt.fun <- cicfamily$wt.fun
	zvec.fun <- cicfamily$zvec.fun
	muhat.fun <- cicfamily$muhat.fun
	ysim.fun <- cicfamily$ysim.fun
	deriv.fun <- cicfamily$deriv.fun
	dev.fun <- cicfamily$dev.fun 
	n <- length(y)
	sm <- 1e-7 
	#sm <- 1e-5
	capl <- length(xmat) / n
	if (capl < 1) {capl <- 0}
	if (round(capl, 8) != round(capl, 1)) {stop ("Incompatible dimensions for xmat!")}
	capk <- length(zmat) / n
	if (capk < 1) {capk <- 0}
	if (round(capk, 8) != round(capk, 1)) {stop ("Incompatible dimensions for zmat!")}
#new:
	capls <- sum(shapes == 17)
####################################################
#get basis functions for the constrained components#
####################################################	
	delta <- NULL
	varlist <- NULL
	xid1 <- NULL; xid2 <- NULL; xpos2 <- 0
	knotsuse <- list(); numknotsuse <- NULL
	mslst <- list()
#new:
	capms <- 0
	if (capl - capls > 0) {
		del1_ans <- makedelta(xmat[, 1], shapes[1], numknots[1], knots[[1]], space = space[1])
		del1 <- del1_ans$amat
		knotsuse[[1]] <- del1_ans$knots 
		mslst[[1]] <- del1_ans$ms
		numknotsuse <- c(numknotsuse, length(del1_ans$knots))
        	m1 <- length(del1) / n
#new code: record the number of columns of del1 if shapes0[1] == 17:
		if (shapes[1] == 17) {capms <- capms + m1}
        	var1 <- 1:m1*0 + 1
		xpos1 <- xpos2 + 1
		xpos2 <- xpos2 + m1
		xid1 <- c(xid1, xpos1)
		xid2 <- c(xid2, xpos2)
		if (capl == 1) {
        		delta <- del1
         		varlist <- var1
          	} else {
	      		for (i in 2:capl) {
#new code:
	        		del2_ans <- makedelta(xmat[,i], shapes[i], numknots[i], knots[[i]], space = space[i])
				del2 <- del2_ans$amat
				knotsuse[[i]] <- del2_ans$knots
				mslst[[i]] <- del2_ans$ms
				numknotsuse <- c(numknotsuse, length(del2_ans$knots))
				m2 <- length(del2) / n
#new code: record the number of columns of del2 if shapes0[i] == 17:
				if (shapes[i] == 17) {capms <- capms + m2}
				xpos1 <- xpos2 + 1
				xpos2 <- xpos2 + m2
				xid1 <- c(xid1, xpos1)
				xid2 <- c(xid2, xpos2)
				delta <- rbind(del1, del2)
				varlist <- 1:(m1 + m2)*0
				varlist[1:m1] <- var1
				varlist[(m1 + 1):(m1 + m2)] <- (1:m2)*0 + i
				var1 <- varlist
				m1 <- m1 + m2
				del1 <- delta
	      		}
	    	}
		if (sum(shapes > 2 & shapes < 5 | shapes > 10 & shapes < 13) > 0 & capk > 0) {
			bigmat <- rbind(1:n*0 + 1, t(zmat), t(xmat[, shapes > 2 & shapes < 5 | shapes > 10 & shapes < 13]), delta)
			np <- 1 + capk + sum(shapes > 2 & shapes < 5 | shapes > 10 & shapes < 13)  + capms
		} else if (sum(shapes > 2 & shapes < 5 | shapes > 10 & shapes < 13) > 0 & capk == 0) {
			bigmat <- rbind(1:n*0 + 1, t(xmat[, shapes > 2 & shapes < 5 | shapes > 10 & shapes < 13]), delta)
			np <- 1 + sum(shapes > 2 & shapes < 5 | shapes > 10 & shapes < 13) + capms
		} else if (sum(shapes > 2 & shapes < 5 | shapes > 10 & shapes < 13) == 0 & capk > 0) {
			bigmat <- rbind(1:n*0 + 1, t(zmat), delta)
			np <- 1 + capk + capms
		} else if (sum(shapes > 2 & shapes < 5 | shapes > 10 & shapes < 13) == 0 & capk == 0) {
			bigmat <- rbind(1:n*0 + 1, delta)
			np <- 1 + capms
		} else {
			print ("error in capk, shapes!")
		} 
#new:
	capm <- length(delta) / n - capms
	} else {
	  	if (capk + capls > 0) {
#new:
			if (capls  < 1 & capk > 0) {
          			bigmat <- rbind(1:n*0 + 1, t(zmat))
          			np <- 1 + capk
			} else if (capls > 0) {
				delta <- NULL; varlist <- NULL
				del1_ans <- makedelta(xmat[,1], 17, numknots[1], knots[[1]], space = space[1])
				del1 <- del1_ans$amat
				knotsuse[[1]] <- del1_ans$knots 
				mslst[[1]] <- del1_ans$ms
				numknotsuse <- c(numknotsuse, length(del1_ans$knots))
				m1 <- length(del1) / n
				var1 <- 1:m1*0 + 1
				xpos1 <- xpos2 + 1
				xpos2 <- xpos2 + m1
				xid1 <- c(xid1, xpos1)
				xid2 <- c(xid2, xpos2)
				if (capls == 1) {
        				delta <- del1
         				varlist <- var1
          			} else {
					for (i in 2:capls) {
	        				del2_ans <- makedelta(xmat[,i], 17, numknots[i], knots[[i]], space = space[i])
						del2 <- del2_ans$amat
						knotsuse[[i]] <- del2_ans$knots
						mslst[[i]] <- del2_ans$ms
						numknotsuse <- c(numknotsuse, length(del2_ans$knots))
						m2 <- length(del2) / n
						xpos1 <- xpos2 + 1
						xpos2 <- xpos2 + m2
						xid1 <- c(xid1, xpos1)
						xid2 <- c(xid2, xpos2)
						delta <- rbind(del1, del2)
						varlist <- 1:(m1 + m2)*0
						varlist[1:m1] <- var1
						varlist[(m1 + 1):(m1 + m2)] <- (1:m2)*0 + i
						var1 <- varlist
						m1 <- m1 + m2
						del1 <- delta
	      				}
				}
				if (capk < 1){
					bigmat <- rbind(1:n*0 + 1, delta)
					capms <- length(delta) / n
					np <- 1 + capms
				} else {
					bigmat <- rbind(1:n*0 + 1, t(zmat), delta)
					capms <- length(delta) / n
					np <- 1 + capk + capms 
				}			
			}
            	} else {bigmat <- matrix(1:n*0 + 1, nrow = 1); capm <- 0; capms <- 0; np <- 1}
	}
	if (!is.null(umbrella.delta)) {
		bigmat <- rbind(bigmat, umbrella.delta)
		capu <- length(umbrella.delta) / n
	} else {capu <- 0}
	if (!is.null(tree.delta)) {
		bigmat <- rbind(bigmat, tree.delta)
		capt <- length(tree.delta) / n
	} else {capt <- 0}
	if (!is.null(umbrella.delta) | !is.null(tree.delta)) 
		delta_ut <- rbind(umbrella.delta, tree.delta)
	if (capl + capk + capu + capt > 0) {
#		if (capl + capu + capt > 0) {
#new:
		if (capl - capls + capu + capt > 0) {
#new:initialize cvec
cvec <- NULL
			if (wt.iter) {
				etahat <- etahat.fun(n, y, fml = family$family)
				gr <- gr.fun(y, etahat, weights, fml = family$family)  
				wt <- wt.fun(etahat, n, weights, fml = family$family)     
				cvec <- wt * etahat - gr
			} else {wt <- wt.fun(etahat, n, weights, fml = family$family)}
			  	zvec <- zvec.fun(cvec, wt, y, fml = family$family)
        		  	gmat <- t(bigmat %*% sqrt(diag(wt)))
			  	dsend <- gmat[, (np + 1):(np + capm + capu + capt), drop = FALSE]
               		  	zsend <- gmat[, 1:np, drop = FALSE] 
				ans <- coneB(zvec, t(dsend), zsend)
			  	etahat <- t(bigmat) %*% ans$coefs
			  	if (wt.iter) {
					muhat <- muhat.fun(etahat, fml = family$family)
			  		diff <- 1
					if (family$family == "binomial") {
						mdiff <- abs(max(muhat) - 1) > sm	
					} else {mdiff <- TRUE}

					nrep <- 0
##########
#iterate!#	
##########
			  		while (diff > sm & mdiff & nrep < n^2){
						oldmu <- muhat	
						nrep <- nrep + 1
						gr <- gr.fun(y, etahat, weights, fml = family$family)	
						wt <- wt.fun(etahat, n, weights, fml = family$family) 
						cvec <- wt * etahat - gr
						#zvec <- cvec / sqrt(wt)
						zvec <- zvec.fun(cvec, wt, y, fml = family$family)						
						gmat <- t(bigmat %*% sqrt(diag(wt)))
						dsend <- gmat[, (np + 1):(np + capm + capu + capt), drop = FALSE]
        	  				zsend <- gmat[, 1:np, drop = FALSE]
						ans <- coneB(zvec, t(dsend), zsend)
						etahat <- t(bigmat) %*% ans$coefs
						muhat <- muhat.fun(etahat, fml = family$family)
						diff <- mean((muhat - oldmu)^2)	
						mdiff <- abs(max(muhat) - 1)
						if (family$family == "binomial") {
							mdiff <- abs(max(muhat) - 1) > sm	
						} else {mdiff <- TRUE}
					}
			 	}
				yhat <- ans$yhat
				coefskeep <- ans$coefs 
########################
#if capk >= 0, we have:#
########################
				zcoefs <- coefskeep[1:(capk + 1)]
######################
#we will always have:#
######################
				vcoefs <- coefskeep[1:np]
#######################
#if capm > 0, we have:#
#######################
				xcoefs <- NULL
				#if (capm > 0) {
				#	xcoefs <- coefskeep[(np + 1):(np + capm)]
				#}
#new:
				if (capl > 0) {
					xcoefs <- coefskeep[(np - capms + 1):(np + capm)]
				}
#######################
#if capu > 0, we have:#
#######################

				ucoefs <- NULL
				if (capu > 0) {
					ucoefs <- coefskeep[(np + 1 + capm):(np + capm + capu)]
				}
#######################
#if capt > 0, we have:#
#######################
				tcoefs <- NULL
				if (capt > 0) {
					tcoefs <- coefskeep[(np + 1 + capm + capu):(np + capm + capu + capt)]
				}

#########################################################		  
#if we have at least one constrained predictor, we have:#
#########################################################

				thvecs <- NULL
				if (capl > 0) {	
#new code:
					dcoefs <- coefskeep[(np - capms + 1):(np + capm)]
					#dcoefs <- coefskeep[(np + 1):(np + capm)]	
	
#####################################################
#thvecs is f(x), where x has one of the eight shapes#
#####################################################

					thvecs <- matrix(nrow = capl, ncol = n)
	    				ncon <- 1
	    				for (i in 1:capl) {
	    	  				thvecs[i,] <- t(delta[varlist == i,]) %*% dcoefs[varlist == i]
#new:
	    	  				#if (shapes[i] > 2 & shapes[i] < 5) { 
						if (shapes[i] > 2 & shapes[i] < 5 | shapes[i] > 10 & shapes[i] < 13) { 
            		    				ncon <- ncon + 1
            		    				#thvecs[i,] <- thvecs[i,] + zcoefs[ncon] * xmat[,i]
							thvecs[i,] <- thvecs[i,] + vcoefs[capk + ncon] * xmat[,i]
            	  				}
	    				}
				}
				thvecs_ut <- NULL
				if (capu + capt > 0) {
					thvecs_ut <- t(delta_ut) %*% coefskeep[(np + 1 + capm):(np + capm + capu + capt)]
				}
				if (!is.null(thvecs_ut)) {
					thvecs <- rbind(thvecs, t(thvecs_ut))
				}
#new:
				llh <-  llh.fun(y, muhat, etahat, n, weights, fml = family$family)
	  			etakeep <- etahat
				muhatkeep <- muhat.fun(etakeep, fml = family$family)
				wtkeep <- wt
#print (sum(abs(coefskeep) > 0))
				df_obs <- sum(abs(coefskeep) > 0)
				if (family$family == "poisson") {		
					mu0 <- mean(y)
					eta0 <- log(mu0)
				} else {mu0 <- NULL}
		} else if (capk + capls > 0 & capl - capls + capt + capu == 0) {
			if (is.null(weights)) {
				weights <- 1:n*0 + 1
			}
			prior.w <- weights
			vmat <- t(bigmat[1:np, , drop = FALSE])			
			if (wt.iter) {
				nrep <- 0
				muhat <- mean(y) + 1:n*0
				etahat <- linkfun(muhat)
				diff <- 1
				if (family$family == "binomial") {
					mdiff <- abs(max(muhat) - 1) > sm	
				} else {mdiff <- TRUE}
				while (diff > sm & mdiff & nrep < n^2) {
					nrep <- nrep + 1
					oldmu <- muhat
					zhat <- etahat + (y - muhat) * deriv.fun(muhat, fml = family$family)				
					#w <- diag(as.vector(prior.w / deriv.fun(muhat)))		
					w <- diag(as.vector(prior.w * (deriv.fun(muhat, fml = family$family))^(-1)))
					b <- solve(t(vmat) %*% w %*% vmat) %*% t(vmat) %*% w %*% zhat
					etahat <- vmat %*% b
					muhat <- muhat.fun(etahat, fml = family$family)		
					diff <- mean((muhat - oldmu)^2)	
					mdiff <- abs(max(muhat) - 1)
					if (family$family == "binomial") {
						mdiff <- abs(max(muhat) - 1) > sm	
					} else {mdiff <- TRUE}
				}
				zcoefs <- b[1:(capk + 1)]
				se.beta <-  sqrt(diag(solve(t(vmat) %*% w %*% vmat)))[1:(capk + 1)]
				zstat <- zcoefs / se.beta
				pvals.beta <-  1 - pchisq(zstat^2, df = 1)
			} else {
				w <- diag(prior.w)
				b <- solve(t(vmat) %*% w %*% vmat) %*% t(vmat) %*% w %*% y
				etahat <- vmat %*% b
				muhat <- muhat.fun(etahat, fml = family$family)	
				sdhat2 <- sum(prior.w * (y - muhat)^2) / (n - np)
				zcoefs <- b[1:(capk + 1)]
				se.beta <-  sqrt(diag(solve(t(vmat) %*% w %*% vmat) * sdhat2))[1:(capk + 1)]
				tstat <- zcoefs / se.beta
				pvals.beta <-  (1 - pt(abs(tstat), df = n - np)) * 2 
			}
#add thvecs if capls > 0:
			thvecs <- NULL
			if (capls > 0) {
				thvecs <- matrix(nrow = capls, ncol = n)
	    			dcoefs <- b[(capk + 2):np]
	    			for (i in 1:capls) {
					thvecs[i,] <- t(delta[varlist == i,]) %*% dcoefs[varlist == i]
	    			}
			}
			llh <-  llh.fun(y, muhat, etahat, n, weights, fml = family$family)
			df_obs <- np
			dfmean <- np
			rslt <- new.env()
			rslt$family <- family 
			rslt$wt.iter <- wt.iter 
			rslt$wt <- diag(w)
			rslt$bigmat <- bigmat
			rslt$etahat <- etahat
			rslt$muhat <- muhat
			rslt$d0 <- np
			rslt$capm <- 0
			rslt$capms <- capms
			rslt$capk <- capk
			rslt$capu <- capu
			rslt$capt <- capt
			rslt$xid1 <- xid1 + np - capms
			rslt$xid2 <- xid2 + np - capms 
			rslt$dfmean <- dfmean	
			rslt$edf0 <- dfmean 
			#rslt$llh <- llh 
			#if (nsim > 0) {
			rslt$cic <- llh + log(1 + 2 * dfmean / (n - np - 1.5 * (dfmean - np)))
			#}
			rslt$zcoefs <- zcoefs
			rslt$coefs <- b
			rslt$vcoefs <- b
			rslt$xcoefs <- b[(capk + 2):np]
			rslt$se.beta <- se.beta 
			rslt$pvals.beta <- pvals.beta 
			rslt$dev <- dev.fun(y, muhat, etahat, weights, fml = family$family)$dev
			rslt$dev.null <- dev.fun(y, muhat, etahat, weights, fml = family$family)$dev.null
			rslt$df <- n - np 
			rslt$df.null <- n - 1
			rslt$resid_df_obs <- n - np - 1.5 * (df_obs - np)
			rslt$vhat <- etahat 		
			rslt$vmat <- vmat	
			rslt$etacomps <- thvecs
			rslt$knots <- knotsuse
			rslt$numknots <- numknotsuse 
			rslt$ms <- mslst
			return (rslt) 
			
		}

##########
#get cic#
##########
		  if (capl - capls + capu + capt > 0 & nsim > 0) {
	  		dfs <- 1:nsim
	  		for (isim in 1:nsim) {
				#set.seed(123)
	    			ysim <- ysim.fun(n, mu0, fml = family$family)
		  		if (wt.iter) {
					etahat <- etahat.fun(n, ysim, fml = family$family)
					gr <- gr.fun(ysim, etahat, weights, fml = family$family)
					wt <- wt.fun(etahat, n, weights, fml = family$family)
					cvec <- wt * etahat - gr
				} else {wt <- wt.fun(etahat, n, weights, fml = family$family)}
					zvec <- zvec.fun(cvec, wt, ysim, fml = family$family)
            				gmat <- t(bigmat %*% sqrt(diag(wt)))
           				dsend <- gmat[, (np + 1):(np + capm + capu + capt), drop = FALSE]
            				zsend <- gmat[, 1:np, drop = FALSE]
					ans <- try(coneB(zvec, t(dsend), zsend))
					if (class(ans) == "try-error") next 
				if (wt.iter) {
						etahat <- t(bigmat) %*% ans$coefs
						muhat <- muhat.fun(etahat, fml = family$family)
						diff <- 1
						if (family$family == "binomial") {
							mdiff <- abs(max(muhat) - 1) > sm	
						} else {mdiff <- TRUE}
##########
#iterate!#
##########
						nrep <- 0 
						while (diff > sm & nrep < n^2 & mdiff > sm) {
							nrep <- nrep + 1
							oldmu <- muhat	
							gr <- gr.fun(ysim, etahat, weights, fml = family$family)
							wt <- wt.fun(etahat, n, weights, fml = family$family)
							cvec <- wt * etahat - gr
							#zvec <- cvec / sqrt(wt)
							zvec <- zvec.fun(cvec, wt, y, fml = family$family)
							gmat <- t(bigmat %*% sqrt(diag(wt)))
							dsend <- gmat[, (np + 1):(np + capm + capu + capt), drop = FALSE]
							zsend <- gmat[, 1:np, drop = FALSE]
							ans <- try(coneB(zvec, t(dsend), zsend))
							if (class(ans) == "try-error") next 
							etahat <- t(bigmat) %*% ans$coefs
							muhat <- muhat.fun(etahat, fml = family$family)
							diff <- mean((muhat - oldmu)^2)	
							if (family$family == "binomial") {
								mdiff <- abs(max(muhat) - 1) > sm	
							} else {mdiff <- TRUE}
						}
				}
	    			dfs[isim] <- sum(abs(ans$coefs) > 0)
	  		}
	  		dfmean <- mean(dfs)
#print (dfmean)
		   } else if (capl - capls + capu + capt > 0 & nsim == 0) {
			dfmean <- NULL
		   } 
###################################################
#if the user does not give any predictor, we have:#
###################################################
	} else {
		rslt <- new.env()
		rslt$family <- family 
		rslt$wt.iter <- wt.iter
		rslt$muhat <- 1:n*0 + mean(y)
		rslt$etahat <- linkfun(muhat)
		rslt$dfmean <- 1
		print ("No predictor is provided")
		return (rslt)	
	}
	if (capl - capls + capu + capt > 0) {
#new:
#exclude the case we only have z or unrestricted smooth
		xid1 <- xid1 + np - capms
		xid2 <- xid2 + np - capms 
		#xid1 <- xid1 + np
		#xid2 <- xid2 + np
		rslt <- new.env()		
		rslt$family <- family 
		rslt$wt.iter <- wt.iter
		rslt$wt <- wtkeep
		rslt$bigmat <- bigmat  
		rslt$etahat <- etakeep
		rslt$muhat <- muhatkeep 
		rslt$d0 <- np
		rslt$capm <- capm
		rslt$capms <- capms 
		rslt$capk <- capk
		rslt$capu <- capu
		rslt$capt <- capt
		rslt$xid1 <- xid1
		rslt$xid2 <- xid2 
		rslt$coefs <- coefskeep
		rslt$vcoefs <- vcoefs
		rslt$xcoefs <- xcoefs
		rslt$zcoefs <- zcoefs
		rslt$ucoefs <- ucoefs
		rslt$tcoefs <- tcoefs 
		rslt$dfmean <- dfmean
		#rslt$llh <- llh
		#if (nsim > 0) {
		if (!is.null(dfmean)) {
#check!
#new: use dfmean - np if n - np - 1.5 * (dfmean - np) < 0
			if ((n - np - 1.5 * (dfmean - np)) <= 0) {
				rslt$cic <- llh + log(1 + 2 * dfmean / (dfmean - np))
			} else {
        			rslt$cic <- llh + log(1 + 2 * dfmean / (n - np - 1.5 * (dfmean - np)))
			}
		} else {rslt$cic <- NULL}
		#if (nsim > 0) {
			#rslt$edf0 <- dfmean - np
		#}
		rslt$edf0 <- dfmean 
		rslt$etacomps <- thvecs
		vmat <- t(bigmat[1:np, , drop = FALSE])
		rslt$vmat <- vmat
		if (is.null(weights)) {
			weights <- 1:n*0 + 1
		}
		prior.w <- weights
		w <- diag(as.vector(prior.w / deriv.fun(muhatkeep, fml = family$family)))
###############################################
#the case capk = 0 and capk >= 0 are combined:#
###############################################
#debugged: vhat -> muhat.fun(vhat)
		vhat <- vmat %*% vcoefs	
		muvhat <- muhat.fun(vhat, fml = family$family)
#debugged: yhat -> muhatkeep
		sse1 <- sum(prior.w * (y - muhatkeep)^2)
		sse0 <- sum(prior.w * (y - muvhat)^2)
#new: use (df_obs - np) if (n - np - 1.5 * (df_obs - np)) <= 0
		if ((n - np - 1.5 * (df_obs - np)) <= 0) {
			sdhat2 <- sse1 / (df_obs - np)
		} else {
			sdhat2 <- sse1 / (n - np - 1.5 * (df_obs - np))
		}
#debugged: vmat -> vmat and duse
#new: coefskeep include zcoefs; bigmat include vmat
		pmat <- vmat
		capbm <- length(bigmat) / n
		bigmat_nv <- bigmat[(np + 1):capbm, , drop = FALSE]	
		coefs_nv <- coefskeep[(np + 1):capbm]	
		duse <- coefs_nv > 1e-8
		if (sum(duse) >= 1) {
			pmat = cbind(vmat, t(bigmat_nv[duse, , drop = FALSE]))
		}
		if (wt.iter) {
			se2 <- solve(t(pmat) %*% w %*% pmat)
		} else {
			se2 <- solve(t(pmat) %*% diag(prior.w) %*% pmat) * sdhat2
		}		 		
		se.beta <- 1:(capk + 1)*0
		tstat <- 1:(capk + 1)*0
		pvals.beta <- 1:(capk + 1)*0
		rslt$zcoefs <- zcoefs 	
		for (i in 1:(capk + 1)) {
			se.beta[i] <- sqrt(se2[i,i])
			tstat[i] <- zcoefs[i] / se.beta[i]
#new code: n - np - 1.5 * (df_obs - np) must be positive
			if ((n - np - 1.5 * (df_obs - np)) <= 0) {
				pvals.beta[i] <- 2 * (1 - pt(abs(tstat[i]),  df_obs - np)) 
				warning ('Effective degrees of freedom is close to the number of observations! Inference about parametric covariates is not reliable!')
			} else {
				pvals.beta[i] <- 2 * (1 - pt(abs(tstat[i]),  n - np - 1.5 * (df_obs - np))) 
			}
		}
		rslt$se.beta <- se.beta
		rslt$pvals.beta <- pvals.beta
		rslt$sse1 <- sse1
		rslt$sse0 <- sse0
		rslt$dev <- dev.fun(y, muhatkeep, etakeep, weights, fml = family$family)$dev
		rslt$dev.null <- dev.fun(y, muhatkeep, etakeep, weights, fml = family$family)$dev.null
		rslt$df <- n - np 
		rslt$df.null <- n - 1
		rslt$resid_df_obs <- n - np - 1.5 * (df_obs - np)
		rslt$df_obs <- df_obs
		rslt$vhat <- vhat
		if (length(knotsuse) == 0) {
			knotsuse <- NULL
		}
		rslt$knots <- knotsuse
		rslt$numknots <- numknotsuse 
		rslt$ms <- mslst
		rslt$capms <- capms
		return (rslt)
	}
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

#######################
#eight shape functions#
####################### 
incr <- function(x, numknots = 0, knots = 0, space = "Q") 
{
    cl <- match.call()
    pars <- match.call()[-1]
    attr(x, "nm") <- deparse(pars$x)
    attr(x, "shape") <- 1
    attr(x, "numknots") <- numknots
    attr(x, "knots") <- knots
    attr(x, "space") <- space
    x
}

decr <- function(x, numknots = 0, knots = 0, space = "Q") 
{
    cl <- match.call()
    pars <- match.call()[-1]
    attr(x, "nm") <- deparse(pars$x)
    attr(x, "shape") <- 2
    attr(x, "numknots") <- numknots
    attr(x, "knots") <- knots
    attr(x, "space") <- space
    x 
} 

conv <- function(x, numknots = 0, knots = 0, space = "Q") 
{
    cl <- match.call()
    pars <- match.call()[-1]
    attr(x, "nm") <- deparse(pars$x)
    attr(x, "shape") <- 3
    attr(x, "numknots") <- numknots
    attr(x, "knots") <- knots
    attr(x, "space") <- space
    x
}

conc <- function(x, numknots = 0, knots = 0, space = "Q") 
{
    cl <- match.call()
    pars <- match.call()[-1]
    attr(x, "nm") <- deparse(pars$x)
    attr(x, "shape") <- 4
    attr(x, "numknots") <- numknots
    attr(x, "knots") <- knots
    attr(x, "space") <- space
    x
}

incr.conv <- function(x, numknots = 0, knots = 0, space = "Q") 
{
    cl <- match.call()
    pars <- match.call()[-1]
    attr(x, "nm") <- deparse(pars$x)
    attr(x, "shape") <- 5
    attr(x, "numknots") <- numknots
    attr(x, "knots") <- knots
    attr(x, "space") <- space
    x
}

decr.conv <- function(x, numknots = 0, knots = 0, space = "Q") 
{
    cl <- match.call()
    pars <- match.call()[-1]
    attr(x, "nm") <- deparse(pars$x)
    attr(x, "shape") <- 6
    attr(x, "numknots") <- numknots
    attr(x, "knots") <- knots
    attr(x, "space") <- space
    x
}

incr.conc <- function(x, numknots = 0, knots = 0, space = "Q") 
{
    cl <- match.call()
    pars <- match.call()[-1]
    attr(x, "nm") <- deparse(pars$x)
    attr(x, "shape") <- 7
    attr(x, "numknots") <- numknots
    attr(x, "knots") <- knots
    attr(x, "space") <- space
    x
}

decr.conc <- function(x, numknots = 0, knots = 0, space = "Q") 
{
    cl <- match.call()
    pars <- match.call()[-1]
    attr(x, "nm") <- deparse(pars$x)
    attr(x, "shape") <- 8
    attr(x, "numknots") <- numknots
    attr(x, "knots") <- knots
    attr(x, "space") <- space
    x
}

s.incr <- function(x, numknots = 0, knots = 0, space = "Q")
{
    cl <- match.call()
    pars <- match.call()[-1]
    attr(x, "nm") <- deparse(pars$x)
    attr(x, "shape") <- 9
    attr(x, "numknots") <- numknots
    attr(x, "knots") <- knots
    attr(x, "space") <- space
    x
}

s.decr <- function(x, numknots = 0, knots = 0, space = "Q") 
{
    cl <- match.call()
    pars <- match.call()[-1]
    attr(x, "nm") <- deparse(pars$x)
    attr(x, "shape") <- 10
    attr(x, "numknots") <- numknots
    attr(x, "knots") <- knots	
    attr(x, "space") <- space
    x
}

s.conv <- function(x, numknots = 0, knots = 0, space = "Q") 
{
    cl <- match.call()
    pars <- match.call()[-1]
    attr(x, "nm") <- deparse(pars$x)
    attr(x, "shape") <- 11
    attr(x, "numknots") <- numknots
    attr(x, "knots") <- knots	
    attr(x, "space") <- space
    x
}

s.conc <- function(x, numknots = 0, knots = 0, space = "Q") 
{
    cl <- match.call()
    pars <- match.call()[-1]
    attr(x, "nm") <- deparse(pars$x)
    attr(x, "shape") <- 12
    attr(x, "numknots") <- numknots
    attr(x, "knots") <- knots	
    attr(x, "space") <- space
    x
}

s.incr.conv <- function(x, numknots = 0, knots = 0, space = "Q") 
{
    cl <- match.call()
    pars <- match.call()[-1]
    attr(x, "nm") <- deparse(pars$x)
    attr(x, "shape") <- 13
    attr(x, "numknots") <- numknots
    attr(x, "knots") <- knots	
    attr(x, "space") <- space
    x
}

s.incr.conc <- function(x, numknots = 0, knots = 0, space = "Q") 
{
    cl <- match.call()
    pars <- match.call()[-1]
    attr(x, "nm") <- deparse(pars$x)
    attr(x, "shape") <- 14
    attr(x, "numknots") <- numknots
    attr(x, "knots") <- knots	
    attr(x, "space") <- space
    x
}

s.decr.conv <- function(x, numknots = 0, knots = 0, space = "Q") 
{
    cl <- match.call()
    pars <- match.call()[-1]
    attr(x, "nm") <- deparse(pars$x)
    attr(x, "shape") <- 15
    attr(x, "numknots") <- numknots
    attr(x, "knots") <- knots	
    attr(x, "space") <- space
    x
}

s.decr.conc <- function(x, numknots = 0, knots = 0, space = "Q") 
{
   cl <- match.call()
   pars <- match.call()[-1]
   attr(x, "nm") <- deparse(pars$x)
   attr(x, "shape") <- 16
   attr(x, "numknots") <- numknots
   attr(x, "knots") <- knots	
   attr(x, "space") <- space
   x
}

s <- function(x, numknots = 0, knots = 0, space = "Q") 
{
   cl <- match.call()
   pars <- match.call()[-1]
   attr(x, "nm") <- deparse(pars$x)
   attr(x, "shape") <- 17
   attr(x, "numknots") <- numknots
   attr(x, "knots") <- knots	
   attr(x, "space") <- space
   x
}

######################################################
#tree function: give the tree shape to x and return x#
######################################################
tree <- function(x)
{
    cl <- match.call()
    pars <- match.call()[-1]
    attr(x, "nm") <- deparse(pars$x)
    attr(x, "shape") <- "tree"
    x

}

##################################################
#tree.fun: make delta to a tree ordering variable#
##################################################
tree.fun <- function(x)
{
    if (min(x) != 0) {
	stop ("A tree ordering variable must have its placebo equal to 0!")
    }
    if (!all(round(x, 0) == x)) {
	stop ("All elements of a tree ordering variable must be integers!")
    }
    #if (any(x < 0))
	#stop ("All elements of a tree ordering variable must be positive!")
    nx <- x
    obs <- 1:length(x)
    delta <- matrix(0, nrow = length(attributes(factor(x))$levels) - 1, ncol = length(x))
    pl <- min(nx)
    for (i in 1:nrow(delta)) {
       nx <- nx[which(nx != pl)]
       pl <- min(nx)
       index <- obs[x == pl]
       delta[i, index] <- 1
    }
  attr(delta, "shape") <- "tree"
  delta
}

##############################################################
#umbrella function: give the umbrella shape to x and return x#
##############################################################
umbrella <- function(x)
{
  cl <- match.call()
  pars <- match.call()[-1]
  attr(x, "nm") <- deparse(pars$x)
  attr(x, "shape") <- "umbrella"
  x
}

###########################################################
#umbrella.fun: make delta to an umbrella ordering variable#
###########################################################
umbrella.fun <- function(x)
{
	amat <- amat.fun(x)
	bmat <- bmat.fun(x)
	constr <- t(rbind(amat, bmat))
	vmat <- qr.Q(qr(constr), complete = TRUE)[, -(1:(qr(constr)$rank)), drop = FALSE]
	if (!is.null(bmat)) {
		wperp <- t(rbind(t(vmat), bmat))
		wmat <- qr.Q(qr(wperp), complete = TRUE)[, -(1:(qr(wperp)$rank)), drop = FALSE]	
		atil <- amat %*% wmat 
		delta <- t(wmat %*% t(atil) %*% solve(atil %*% t(atil)))
	} else {
		delta <- t(t(amat) %*% solve(amat %*% t(amat)))
	}
	attr(delta, "shape") <- "umbrella"
	delta
}

###################################################
#find delta for a specific predictor x and a shape#
################################################### 
makedelta = function(x, sh, numknots = 0, knots = 0, space = "E", suppre = FALSE, interp = FALSE) {
	n = length(x)
# find unique x values
#round(x,8) will make 0 edge in amat!
	#xu = sort(unique(round(x, 8)))
	xu = sort(unique(x))
	n1 = length(xu)
	sm = 1e-7
	ms = NULL
#  increasing or decreasing
	if (sh < 3) {
		amat = matrix(0, nrow = n1 - 1, ncol = n)
		for (i in 1: (n1 - 1)) {
			amat[i, x > xu[i]] = 1
		}
		if (sh == 2) {amat = -amat}
		for (i in 1:(n1 - 1)) {amat[i, ] = amat[i, ] - mean(amat[i, ])}
	} else if (sh == 3 | sh == 4) {
#  convex or concave
		amat = matrix(0, nrow = n1 - 2 ,ncol = n)
		for (i in 1: (n1 - 2)) {
			amat[i, x > xu[i]] = x[x > xu[i]] - xu[i]
		}
		if (sh == 4) {amat = -amat}
		xm = cbind(1:n*0+1,x)
		xpx = solve(t(xm) %*% xm)
		pm = xm %*% xpx %*% t(xm)
		dmat = amat - amat %*% t(pm)
	} else if (sh > 4 & sh < 9) {
		amat = matrix(0, nrow = n1 - 1, ncol = n)
		if (sh == 5) { ### increasing convex
			for (i in 1:(n1 - 1)) {
				amat[i, x > xu[i]] = (x[x > xu[i]] - xu[i]) / (max(x) - xu[i])
			}
			for (i in 1:(n1 - 1)) {amat[i,] = amat[i,] - mean(amat[i,])}
		} else if (sh == 6) {  ## decreasing convex
			for (i in 1:(n1 - 1)) {
				amat[i, x < xu[i + 1]] = (x[x < xu[i + 1]] - xu[i + 1]) / (min(x) - xu[i + 1])
			}
			for (i in 1:(n1 - 1)) {amat[i,] = amat[i,] - mean(amat[i,])}
		} else if (sh == 7) { ## increasing concave
			for (i in 1:(n1 - 1)) {
				amat[i, x < xu[i + 1]] = (x[x < xu[i + 1]] - xu[i + 1]) / (min(x) - xu[i + 1])
			}
			for (i in 1:(n1 - 1)) {amat[i,] = -amat[i,] + mean(amat[i,])}		
		} else if (sh == 8) {## decreasing concave
			for (i in 1:(n1 - 1)) {
				amat[i, x > xu[i]] = (x[x > xu[i]] - xu[i]) / (max(x) - xu[i])
			}
			for (i in 1:(n1 - 1)) {amat[i,] = -amat[i,] + mean(amat[i,])}
		}
	} else if (sh > 8 & sh < 18) {
		#if (all(knots == 0) & numknots == 0) {
		if (length(knots) < 2 & numknots == 0) {
			if (sh == 9 | sh == 10) {#1 2
				k = trunc(n1^(1/5)) + 4
			} else {k = trunc(n1^(1/7) + 4)}
			if (space == "Q") {
				t = quantile(xu, probs = seq(0, 1, length = k), names = FALSE)
			}
			if (space == "E") {
				t = 0:k / k * (max(x) - min(x)) + min(x)
			} 
		#} else if (any(knots != 0) & numknots == 0) {
		} else if (length(knots) >= 2 & numknots == 0) {
			t = knots
		#} else if (all(knots == 0) & numknots != 0) {
		} else if (length(knots) < 2 & numknots != 0) {
			if (space == "Q") {
				t = quantile(xu, probs = seq(0, 1, length = numknots), names = FALSE)
			} 
			if (space == "E") {
				k = numknots
				#if (sh == 9 | sh == 10) {#1 2
				#	k = trunc(n1^(1/5)) + 4
				#} else {k = trunc(n1^(1/7) + 4)}
				t = 0:k / k * (max(x) - min(x)) + min(x)
			}
		#} else if (any(knots != 0) & numknots != 0) {
		} else if (length(knots) >= 2 & numknots != 0) {
			#t0 = quantile(xu, probs = seq(0, 1, length = numknots), names = FALSE)
			t = knots
			if (!suppre) {
				print("'knots' is used! 'numknots' is not used!")
			}
			#print ("'knots' is used!")
			#if (numknots != length(knots)) {
			#	if (!suppre) {
			#		print("length(knots) is not equal to 'numknots'! 'knots' is used!")
			#	}
			#} else if (any(t0 != knots)) {
			#	if (!suppre) {
			#		print("equal x-quantiles knots != 'knots'! 'knots' is used! ") 
			#	}
			#}
		}
		if (sh == 9) {#1			
			amat_ans = monincr(x, t, interp)
			amat = amat_ans$sigma
			ms = amat_ans$ms
		} else if (sh == 10) {#2
			amat_ans = mondecr(x, t, interp)
			amat = amat_ans$sigma
			ms = amat_ans$ms
		} else if (sh == 11) {#3
			amat_ans = convex(x, t, interp)
			amat = amat_ans$sigma
			ms = amat_ans$ms
		} else if (sh == 12) {#4
			amat_ans = concave(x, t, interp)
			amat = amat_ans$sigma
			ms = amat_ans$ms
		} else if (sh == 13) {#5
			amat_ans = incconvex(x, t, interp)
			amat = amat_ans$sigma
			ms = amat_ans$ms
		} else if (sh == 14) {#6
			amat_ans = incconcave(x, t, interp)
			amat = amat_ans$sigma
			ms = amat_ans$ms
		} else if (sh == 15) {#7
			#amat_ans = -incconcave(x, t, interp)
			amat_ans = incconcave(x, t, interp)
			amat = -amat_ans$sigma
			if (!interp) {
				ms = -amat_ans$ms
			}
		} else if (sh == 16) {#8
			#amat_ans = -incconvex(x, t, interp)
			amat_ans = incconvex(x, t, interp)
			amat = -amat_ans$sigma
			if (!interp) {
				ms = -amat_ans$ms
			}
		} else if (sh == 17) {#unconstrained
			amat_ans = incconvex(x, t, interp)
			amat = amat_ans$sigma
			ms = amat_ans$ms
			#amat = -incconcave(x, t)
			#amat = rbind(x, t(bcspl(x, m = length(t), knots = t)$bmat)) 
			#amat = rbind(x, convex(x, t))
		}
	}
	#if (sh < 9) {
	#	rslt = list(amat = amat, knots = 0, ms = ms)
	#} else {
	#	rslt = list(amat = amat, knots = t, ms = ms)
	#}
	if (sh < 9) {t = 0}	
	rslt = list(amat = amat, knots = t, ms = ms)
	rslt
}

# Monotone increasing
monincr = function(xs, t, interp = FALSE) {
	n = length(xs)
	x = sort(xs)
	k = length(t) - 2
	m = k + 2
	sigma = matrix(0, nrow = m, ncol = n)
	obs = 1:n
	knt = 1:m
	for (i in 1:(k+2)) {knt[i] = min(obs[abs(x - t[i]) == min(abs(x - t[i]))])}
	for (j in 1:(k-1)) {
		index = x >= t[1] & x <= t[j]
		sigma[j, index] = 0

		index = x > t[j] & x <= t[j+1]
		sigma[j, index] = (x[index] - t[j])^2 / (t[j+2] - t[j]) / (t[j+1] - t[j])

		index = x > t[j+1] & x <= t[j+2]
		sigma[j, index] = 1 - (x[index] - t[j+2])^2 / (t[j+2] - t[j+1]) / (t[j+2] - t[j])
	    
		index = x > t[j+2] #& x <= t[m]
		sigma[j, index] = 1
	}
	index = x >= t[1] & x <= t[k]
	sigma[k, index] = 0
	
	index = x > t[k] & x <= t[k+1]
	sigma[k, index] = (x[index] - t[k])^2 / (t[k+2] - t[k]) / (t[k+1] - t[k])
	
	index = x > t[k+1] & x <= t[k+2]
	sigma[k, index] = 1 - (x[index] - t[k+2])^2 / (t[k+2] - t[k+1]) / (t[k+2] - t[k])
	
	index = x >= t[1] & x <= t[2]
	sigma[k+1, index] = 1 - (t[2] - x[index])^2 / (t[2] - t[1])^2

	index = x > t[2] 
	sigma[k+1, index] = 1
	
	index = x >= t[1] & x <= t[k+1]
	sigma[k+2, index] = 0
	
	index = x > t[k+1] & x <= t[k+2]
	sigma[k+2, index] = (x[index] - t[k+1])^2 / (t[k+2] - t[k+1])^2
	
#new:
	ms = NULL
	if (!interp) {
		ms = apply(sigma, 1, mean)
		for (i in 1:m) {
			sigma[i,] = sigma[i,] - mean(sigma[i,])
			sigma[i,] = sigma[i, rank(xs)]
		} 
	} else {
		for (i in 1:m) {
			#sigma[i,] = sigma[i,] - mean(sigma[i,])
			sigma[i,] = sigma[i, rank(xs)]
		} 
	}
	rslt = list(sigma = sigma, ms = ms)
	rslt
}	


########################################################
# Monotone decreasing
mondecr = function(xs, t, interp = FALSE) {
	x = sort(xs)
	n = length(x)
	k = length(t) - 2
	m = k + 2
	sigma = matrix(0, nrow = m, ncol = n)
	obs = 1:n
	#knt = 1:m
	#for (i in 1:(k + 2)) {knt[i] = min(obs[abs(x - t[i]) == min(abs(x - t[i]))])}
	#t = x[knt]
	for (j in 1:(k - 1)) {
	 	index = x >= t[1] & x <= t[j]
	 	sigma[j, index] = 1

		index = x > t[j] & x <= t[j+1]
	 	sigma[j, index] = 1 - (x[index] - t[j])^2 / (t[j+2] - t[j]) / (t[j+1] - t[j])

	    	index = x > t[j+1] & x <= t[j+2]
	    	sigma[j, index] = (x[index] - t[j+2])^2 / (t[j+2] - t[j+1]) / (t[j+2] - t[j])

	    	index = x > t[j+2] 
	    	sigma[j, index] = 0
	}

	index = x >= t[1] & x <= t[k]
	sigma[k, index] = 1
	
	index = x > t[k] & x <= t[k+1]
	sigma[k, index] = 1 - (x[index] - t[k])^2 / (t[k+2] - t[k]) / (t[k+1] - t[k])

	index = x > t[k+1] & x <= t[k+2]
	sigma[k, index] = (x[index] - t[k+2])^2 / (t[k+2] - t[k+1]) / (t[k+2] - t[k])

	index = x >= t[1] & x <= t[2]
	sigma[k+1, index] = (t[2] - x[index])^2 / (t[2] - t[1])^2

	index = x > t[2] 
	sigma[k+1, index] = 0

	index = x >= t[1] & x <= t[k+1]
	sigma[k+2, index] = 1
	
	index = x > t[k+1] & x <= t[k+2]
	sigma[k+2, index] = 1 - (x[index] - t[k+1])^2 / (t[k+2] - t[k+1])^2

	ms = NULL
	if (!interp) {
		ms = apply(sigma, 1, mean)
		for (i in 1:m) {
			sigma[i,] = sigma[i,] - mean(sigma[i,])
			sigma[i,] = sigma[i, rank(xs)]
		} 
	} else {
		for (i in 1:m) {
			#sigma[i,] = sigma[i,] - mean(sigma[i,])
			sigma[i,] = sigma[i, rank(xs)]
		} 
	}
	rslt = list(sigma = sigma, ms = ms)
	rslt
}

########################################################
# Convex
convex = function(xs, t, interp = FALSE) {
	x = sort(xs)
	n = length(x)
	k = length(t) - 2
	m = k + 2
	sigma = matrix(0, nrow = m, ncol = n)
	obs = 1:n
	#knt = 1:m
	#for (i in 1:(k+2)) {knt[i] = min(obs[abs(x - t[i]) == min(abs(x - t[i]))])}
	for (j in 1:(k-1)) {
	 	index = x >= t[1] & x <= t[j]
	 	sigma[j, index] = 0
	 	
	 	index = x > t[j] & x <= t[j+1]
	 	sigma[j, index] = (x[index] - t[j])^3 / (t[j+2] - t[j]) / (t[j+1] - t[j]) / 3
	    
	   	index = x > t[j+1] & x <= t[j+2]
	    	sigma[j, index] = x[index] - t[j+1] - (x[index] - t[j+2])^3 / (t[j+2] - t[j]) / (t[j+2] - t[j+1]) / 3 + (t[j+1] - t[j])^2 / 3 /(t[j+2] - t[j]) - (t[j+2] - t[j+1])^2 / 3 / (t[j+2] - t[j])

	    	index = x > t[j+2]
	    	sigma[j, index] = (x[index] - t[j+1]) + (t[j+1] - t[j])^2 / 3 / (t[j+2] - t[j]) - (t[j+2] - t[j+1])^2 / 3 / (t[j+2] - t[j])
	}
	index = x >= t[1] & x <= t[k]
	sigma[k, index] = 0
	
	index = x > t[k] & x <= t[k+1]
	sigma[k, index] = (x[index] - t[k])^3 / (t[k+2] - t[k]) / (t[k+1] - t[k]) / 3

	index = x > t[k+1] & x <= t[k+2]
	sigma[k, index] = x[index] - t[k+1] - (x[index] - t[k+2])^3 / (t[k+2] - t[k]) / (t[k+2] - t[k+1]) / 3 + (t[k+1] - t[k])^2 / 3 / (t[k+2] -t[k]) - (t[k+2] - t[k+1])^2 / 3 / (t[k+2] - t[k])
	
	index = x >= t[1] & x <= t[2]
	sigma[k+1, index] = x[index] - t[1] + (t[2] - x[index])^3 / (t[2] - t[1])^2 / 3 - (t[2] - t[1]) / 3
	
	index = x > t[2] 
	sigma[k+1, index] = x[index] - t[1] - (t[2] - t[1]) / 3
		
	index = x >= t[1] & x <= t[k+1]
	sigma[k+2, index] = 0
	
	index = x > t[k+1] & x <= t[k+2]
	sigma[k+2, index] = (x[index] - t[k+1])^3 / (t[k+2] - t[k+1])^2 / 3
	
	ms = NULL
	if (!interp) {
		xm = cbind(1:n*0+1, x)
		pm = xm %*% solve(t(xm) %*% xm) %*% t(xm)
		ms = matrix(0, nrow = nrow(sigma), ncol = ncol(sigma))
		for (i in 1:m) {
			ms[i,] = pm %*% sigma[i,]
			ms[i,] = ms[i, rank(xs)]		
			sigma[i,] = sigma[i,] - pm %*% sigma[i,]
			sigma[i,] = sigma[i, rank(xs)]
		}
	} else {
		for (i in 1:m) {
			#sigma[i,] = sigma[i,] - pm %*% sigma[i,]
			sigma[i,] = sigma[i, rank(xs)]
		}
	}
	rslt = list(sigma = sigma, ms = ms)
	rslt
}


########################################################
# Concave
concave = function(xs, t, interp = FALSE) {
	x = sort(xs)
	n = length(x)
	k = length(t) - 2
	m = k + 2
	sigma = matrix(0, nrow = m, ncol = n)
	obs = 1:n
	#knt = 1:m
	#for (i in 1:(k+2)) {knt[i] = min(obs[abs(x - t[i]) == min(abs(x - t[i]))])}
	#t = x[knt]
	for (j in 1:k) {
	 	index = x >= t[1] & x <= t[j]
	 	sigma[j, index] = x[index] - t[1]
	 	
	 	index = x > t[j] & x <= t[j+1]
	 	sigma[j, index] = t[j] - t[1] + ((t[j+1] - t[j])^3 - (t[j+1] - x[index])^3) / 3 / (t[j+1] - t[j]) / (t[j+2] - t[j]) + (x[index] - t[j]) * (t[j+2] - t[j+1]) / (t[j+2] - t[j])
	    
	        index = x > t[j+1] & x <= t[j+2]
	    	sigma[j, index] = t[j] - t[1] + (t[j+1] - t[j])^2 / 3 / (t[j+2] - t[j]) + (t[j+2] - t[j+1]) * (t[j+1] - t[j]) / (t[j+2] - t[j]) + ((t[j+2] - t[j+1])^3 - (t[j+2] - x[index])^3) / 3 / (t[j+2] - t[j+1]) / (t[j+2] - t[j])	
 	   
 	   	index = x > t[j+2]
 	   	sigma[j, index] = t[j] - t[1] + (t[j+1] - t[j])^2 / 3 / (t[j+2] - t[j]) + (t[j+2] - t[j+1]) * (t[j+1] - t[j]) / (t[j+2] - t[j]) + (t[j+2] - t[j+1])^2 / 3 / (t[j+2] - t[j])
	}

	index = x >= t[1] & x <= t[2]
	sigma[k+1, index] = -(t[2] - x[index])^3 / 3 / (t[2] - t[1])^2
	
	index = x > t[2] 
	sigma[k+1, index] = 0
	
	index = x >= t[1] & x <= t[k+1]
	sigma[k+2, index] = x[index] - t[1]
	
	index = x > t[k+1] & x <= t[k+2]
	sigma[k+2, index] = t[k+1] - t[1] + ((t[k+2] - t[k+1])^2 * (x[index] - t[k+1]) - (x[index] - t[k+1])^3 / 3) / (t[k+2] - t[k+1])^2
	
	ms = NULL
	if (!interp) {
		xm = cbind(1:n*0+1, x)
		pm = xm %*% solve(t(xm) %*% xm) %*% t(xm)
		ms = matrix(0, nrow = nrow(sigma), ncol = ncol(sigma))
		for (i in 1:m) {
			ms[i,] = pm %*% sigma[i,]
			ms[i,] = ms[i, rank(xs)]		
			sigma[i,] = sigma[i,] - pm %*% sigma[i,]
			sigma[i,] = sigma[i, rank(xs)]
		}
	} else {
		for (i in 1:m) {
			#sigma[i,] = sigma[i,] - pm %*% sigma[i,]
			sigma[i,] = sigma[i, rank(xs)]
		}
	}	
	rslt = list(sigma = sigma, ms = ms)
	rslt
}

########################################################
# Increasing and Convex
incconvex = function(xs, t, interp = FALSE) {
	x = sort(xs)
	n = length(x)
	k = length(t) - 2
	m = k + 3
	sigma = matrix(0, nrow = m, ncol = n)
	obs = 1:n
	#knt = 1:(k+2)
	#for (i in 1:(k+2)) {knt[i] = min(obs[abs(x - t[i]) == min(abs(x - t[i]))])}
	for (j in 1:(k-1)) {
	 	index = x >= t[1] & x <= t[j]
	 	sigma[j, index] = 0
	 	
	 	index = x > t[j] & x <= t[j+1]
	 	sigma[j, index] = (x[index] - t[j])^3 / (t[j+2] - t[j]) / (t[j+1] - t[j]) / 3
	    
	    	index = x > t[j+1] & x <= t[j+2]
	    	sigma[j, index] = x[index] - t[j+1] - (x[index] - t[j+2])^3 / (t[j+2] - t[j]) / (t[j+2] - t[j+1]) / 3 + (t[j+1] - t[j])^2 / 3 /(t[j+2] - t[j]) - (t[j+2] - t[j+1])^2 / 3 / (t[j+2] - t[j])
	    
	  	index = x > t[j+2] 
	    	sigma[j, index] = (x[index] - t[j+1]) + (t[j+1] - t[j])^2 / 3 / (t[j+2] - t[j]) - (t[j+2] - t[j+1])^2 / 3 / (t[j+2] - t[j])
	}
	index = x >= t[1] & x <= t[k]
	sigma[k, index] = 0
	
	index = x > t[k] & x <= t[k+1]
	sigma[k, index] = (x[index] - t[k])^3 / (t[k+2] - t[k]) / (t[k+1] - t[k]) / 3
	   
	index = x > t[k+1] & x <= t[k+2]
	sigma[k, index] = x[index] - t[k+1] - (x[index] - t[k+2])^3 / (t[k+2] - t[k]) / (t[k+2] - t[k+1]) / 3 + (t[k+1] - t[k])^2 / 3 / (t[k+2] - t[k]) - (t[k+2] - t[k+1])^2 / 3 / (t[k+2] - t[k])
	
	index = x >= t[1] & x <= t[2]
	sigma[k+1, index] = x[index] - t[1] + (t[2] - x[index])^3 / (t[2] - t[1])^2 / 3 - (t[2] - t[1]) / 3
	
	index = x > t[2] 
	sigma[k+1, index] = x[index] - t[1] - (t[2] - t[1]) / 3
	
	index = x >= t[1] & x <= t[k+1]
	sigma[k+2, index] = 0
	
	index = x > t[k+1] & x <= t[k+2]
	sigma[k+2, index] = (x[index] - t[k+1])^3 / (t[k+2] - t[k+1])^2 / 3
	
	sigma[k+3,] = x

	ms = NULL
	if (!interp) {
		ms = apply(sigma, 1, mean)
		for (i in 1:m) {
			sigma[i,] = sigma[i,] - mean(sigma[i,])
			sigma[i,] = sigma[i, rank(xs)]
		} 
	} else {
		for (i in 1:m) {
			#sigma[i,] = sigma[i,] - mean(sigma[i,])
			sigma[i,] = sigma[i, rank(xs)]
		} 
	}
	rslt = list(sigma = sigma, ms = ms)
	rslt
}

########################################################
# Increasing and Concave
incconcave = function(xs, t, interp = FALSE) {
	x = sort(xs)
	n = length(x)
	k = length(t) - 2
	m = k + 3
	sigma = matrix(0, nrow = m, ncol = n)
	obs = 1:n
	#knt = 1:(k+2)
	#for(i in 1:(k+2)) {knt[i] = min(obs[abs(x - t[i]) == min(abs(x - t[i]))])}
	for (j in 1:k) {
	 	index = x >= t[1] & x <= t[j]
	 	sigma[j, index] = x[index] - t[1]
	 	
	 	index = x > t[j] & x <= t[j+1]
	 	sigma[j, index] = t[j] - t[1] + ((t[j+1] - t[j])^3 - (t[j+1] - x[index])^3) / 3 / (t[j+1] - t[j]) / (t[j+2] - t[j]) + (x[index] - t[j]) * (t[j+2] - t[j+1]) / (t[j+2] - t[j])
	    
	    	index = x > t[j+1] & x <= t[j+2]
	    	sigma[j, index] = t[j] - t[1] + (t[j+1] - t[j])^2 / 3 / (t[j+2] - t[j]) + (t[j+2] - t[j+1]) * (t[j+1] - t[j]) / (t[j+2] - t[j]) + ((t[j+2] - t[j+1])^3 - (t[j+2] - x[index])^3) / 3 / (t[j+2] - t[j+1]) / (t[j+2] - t[j])
	   
	    	index = x > t[j+2] 
	    	sigma[j, index] = t[j] - t[1] + (t[j+1] - t[j])^2 / 3 / (t[j+2] - t[j]) + (t[j+2] - t[j+1]) * (t[j+1] - t[j]) / (t[j+2] - t[j]) + (t[j+2] - t[j+1])^2 / 3 /(t[j+2] - t[j])
	    
	}

	index = x >= t[1] & x <= t[2]
	sigma[k+1, index] = -(t[2] - x[index])^3 / 3 / (t[2] - t[1])^2
	
	index = x > t[2] 
	sigma[k+1, index] = 0
	
	index = x >= t[1] & x <= t[k+1]
	sigma[k+2, index] = x[index] - t[1]
	
	index = x > t[k+1] & x <= t[k+2]
	sigma[k+2, index] = t[k+1] - t[1] + ((t[k+2] - t[k+1])^2 * (x[index] - t[k+1]) - (x[index] - t[k+1])^3 / 3) / (t[k+2] - t[k+1])^2
	
	sigma[k+3, ] = x

	ms = NULL
	if (!interp) {
		ms = apply(sigma, 1, mean)
		for (i in 1:m) {
			sigma[i,] = sigma[i,] - mean(sigma[i,])
			sigma[i,] = sigma[i, rank(xs)]
		} 
	} else {
		for (i in 1:m) {
			#sigma[i,] = sigma[i,] - mean(sigma[i,])
			sigma[i,] = sigma[i, rank(xs)]
		} 
	}
	rslt = list(sigma = sigma, ms = ms)
	rslt
}


##############
#summary.cgam#
##############
summary.cgam <- function(object,...) {
	if (!is.null(object$zcoefs)) {
		family <- object$family 
		resid_df_obs <- object$resid_df_obs
		wt.iter <- object$wt.iter
		coefs <- object$zcoefs
		se <- object$se.beta
		tval <- coefs / se
		pvalbeta <- object$pvals.beta
		n <- length(coefs)
		sse0 <- object$SSE0
		sse1 <- object$SSE1
		cic <- object$cic
		deviance <- object$deviance
		null_deviance <- object$null_deviance
		df <- object$df
		null_df <- object$null_df
		zid <- object$zid
		#zid1 <- object$zid1 - 1 - length(shapes)
		#zid2 <- object$zid2 - 1 - length(shapes)
#new: zid1, zid2 just index zmat not bigmat
		zid1 <- object$zid1
		zid2 <- object$zid2
		tms <- object$tms
		is_param <- object$is_param
		is_fac <- object$is_fac
		vals <- object$vals
		if (wt.iter) {
			rslt1 <- data.frame("Estimate" = round(coefs, 4), "StdErr" = round(se, 4), "z.value" = round(tval, 4), "p.value" = round(pvalbeta, 4))			
			rownames(rslt1)[1] <- "(Intercept)"
			if (n > 1) {
				lzid <- length(zid1)
				for (i in 1:lzid) {
					pos1 <- zid1[i]; pos2 <- zid2[i]
					for (j in pos1:pos2) {
						if (!is_param[i]) {
							rownames(rslt1)[j + 1] <- paste(attributes(tms)$term.labels[zid[i] - 1], rownames(rslt1)[j + 1], sep = "")						
						} else {
							rownames(rslt1)[j + 1] <- paste(attributes(tms)$term.labels[zid[i] - 1], vals[j], sep = "")						
						}	
					}
				}
			}
			rslt1 <- as.matrix(rslt1)
		} else {
			rslt1 <- data.frame("Estimate" = round(coefs, 4), "StdErr" = round(se, 4), "t.value" = round(tval, 4), "p.value" = round(pvalbeta, 4))
			rownames(rslt1)[1] <- "(Intercept)"
			if (n > 1) {
				lzid <- length(zid1)
				for (i in 1:lzid) {
					pos1 <- zid1[i]; pos2 <- zid2[i]; 
					for (j in pos1:pos2) {
						if (!is_param[i]) {
							rownames(rslt1)[j + 1] <- paste(attributes(tms)$term.labels[zid[i] - 1], rownames(rslt1)[j + 1], sep = "")						
						} else {
							rownames(rslt1)[j + 1] <- paste(attributes(tms)$term.labels[zid[i] - 1], vals[j], sep = "")						
						}	
					}
				}
			}
			rslt1 <- as.matrix(rslt1)
		}
		#if (!is.null(sse0) & !is.null(sse1)) {
			#rslt2 <- cbind(SSE.Linear = sse0, SSE.Full = sse1)
#new:
		#	rslt2 <- data.frame("SSE.Linear" = sse0, "SSE.Full" = sse1)
		#	rownames(rslt2)[1] <- ""
		#	ans <- list(call = object$call, coefficients = rslt1, residuals = rslt2, zcoefs = coefs, cic = cic, null_deviance = null_deviance, null_df = null_df, deviance = deviance, df = df, resid_df_obs = resid_df_obs, family = family) 
		#	class(ans) <- "summary.cgam"
		#	ans
		#} else {
			ans <- list(call = object$call, coefficients = rslt1, zcoefs = coefs, cic = cic, null_deviance = null_deviance, null_df = null_df, deviance = deviance, df = df, resid_df_obs = resid_df_obs, family = family)
			class(ans) <- "summary.cgam"
			ans
		#}
	} else {
		ans <- list(zcoefs = object$zcoefs)
		class(ans) <- "summary.cgam"
		ans
	}
}

#############
#summary.wps#
#############
summary.wps <- function(object,...) {
	if (!is.null(object$zcoefs)) {
		coefs <- object$zcoefs
		se <- object$se.beta
		#tval <- object$tz
		pvalbeta <- object$pvals.beta
		tval <- coefs / se
		n <- length(coefs)
		#sse0 <- object$SSE0
		#sse1 <- object$SSE1
		zid <- object$zid
#new: zid1, zid2 just index zmat not bigmat
		zid1 <- object$zid1
		zid2 <- object$zid2
		#tms <- object$tms
		#zmat <- object$zmat
		#is_mat <- object$is_mat
		is_param <- object$is_param
		is_fac <- object$is_fac
		vals <- object$vals
		tms <- object$tms
		rslt1 <- data.frame("Estimate" = round(coefs, 4), "StdErr" = round(se, 4), "t.value" = round(tval, 4), "p.value" = round(pvalbeta, 4))
		rownames(rslt1)[1] <- "(Intercept)"
		if (n > 1) {
			lzid <- length(zid1)
			for (i in 1:lzid) {
				pos1 <- zid1[i]; pos2 <- zid2[i]
				for (j in pos1:pos2) {
					if (!is_param[i]) {
						rownames(rslt1)[j + 1] <- paste(attributes(tms)$term.labels[zid[i] - 1], rownames(rslt1)[j + 1], sep = "")						
					} else {
						rownames(rslt1)[j + 1] <- paste(attributes(tms)$term.labels[zid[i] - 1], vals[j], sep = "")					
					}	
				}
			}
		}
		rslt1 <- as.matrix(rslt1)
		#if (!is.null(sse0) & !is.null(sse1)) {
		#	rslt2 <- data.frame("SSE.Linear" = sse0, "SSE.Full" = sse1)
#new:
		#	rownames(rslt2)[1] <- ""
			#rslt2 <- as.matrix(rslt2)
		#	ans <- list(call = object$call, coefficients = rslt1, residuals = rslt2, zcoefs = coefs) 
		#	class(ans) <- "summary.wps"
		#	ans
		#} else {
			ans <- list(call = object$call, coefficients = rslt1, zcoefs = coefs)
			class(ans) <- "summary.wps"
			ans
		#}
	} else {
		ans <- list(zcoefs = object$zcoefs)
		class(ans) <- "summary.wps"
		ans
	}
}


#####################
#print.summary.cgam #
#####################
print.summary.cgam <- function(x,...) {
	if (!is.null(x$zcoefs)) {
	#if (!is.null(x$se.beta)) {
		cat("Call:\n")
		print(x$call)
		cat("\n")
		cat("Coefficients:")
		cat("\n")
		printCoefmat(x$coefficients, P.values = TRUE, has.Pvalue = TRUE)
		cat("\n")
		if (x$family$family == "binomial") {
			message("(Dispersion parameter for binomial family taken to be 1)", "\n")
		}
		if (x$family$family == "poisson") {
			message("(Dispersion parameter for poisson family taken to be 1)", "\n")
		} 
		if (x$family$family == "gaussian") {
			message("(Dispersion parameter for gaussian family taken to be ", round(x$deviance/x$df,4), " )", "\n")
		} 
		message("Null deviance: ", round(x$null_deviance,4), " ", "on ", x$null_df, " ", "degrees of freedom")
		message("Residual deviance: ", round(x$deviance,4), " ", "on ", x$resid_df_obs, " ", "observed degrees of freedom")
		#if (is.null(x$cic)) {
		#	message("CIC value is not available when there is no shape-restricted predictor")
		#} else {message("CIC: ", round(x$cic,4))}
		if (!is.null(x$cic)) {
			message("CIC: ", round(x$cic,4))
		}
		#if (!is.null(x$residuals)) {
		#	cat("==============================================================", "\n")
		#	cat("Call:\n")
		#	print(x$call)
		#	cat("\n")
		#	printCoefmat(x$residuals, P.values = TRUE, has.Pvalue = TRUE)
		#}
	} else {
		print ("No linear predictor is defined")	
		#print ("Residual degree of freedom is negative!")
	}
}

####################
#print.summary.wps #
####################
print.summary.wps <- function(x,...) {
	if (!is.null(x$zcoefs)) {
	#if (!is.null(x$se.beta)) {
		cat("Call:\n")
		print(x$call)
		cat("\n")
		cat("Coefficients:")
		cat("\n")
		printCoefmat(x$coefficients, P.values = TRUE, has.Pvalue = TRUE)
		#if (!is.null(x$residuals)) {
			#cat("==============================================================", "\n")
		#	cat("+--------------------------------------+\n")
  		#	cat("|         SSE.Linear vs SSE.Full       |\n")
#printCoefmat(x$residuals, P.values = TRUE, has.Pvalue = TRUE)
  		#	cat("+--------------------------------------+\n")
			#cat("Call:\n")
			#print(x$call)
			#cat("\n")
		#	printCoefmat(x$residuals, P.values = TRUE, has.Pvalue = TRUE)
		#}
	} else {
		print ("No linear predictor is defined")	
	}
}


##############
#predict.cgam#
##############
predict.cgam = function(object, newData,...) {
	family = object$family
	cicfamily = CicFamily(family)
	muhat.fun = cicfamily$muhat.fun	
	shapes = object$shapes
	np = object$d0; capm = object$capm; capk = object$capk; capt = object$capt; capu = object$capu
#new:
	xid10 = object$xid1; xid20 = object$xid2; 
	uid1 = object$uid1; uid2 = object$uid2; tid1 = object$tid1; tid2 = object$tid2 
#new:
	xmat0 = object$xmat; knots0 = object$knots; numknots0 = object$numknots; sps0 = object$sps; ms0 = object$ms; xmatnms = object$xmatnms; capms = object$capms
	zmat = object$zmat; umb = object$umb; tr = object$tr
	bigmat = object$bigmat; umbrella.delta = object$umbrella.delta; tree.delta = object$tree.delta
	coefs = object$coefs; zcoefs = object$zcoefs; vcoefs = object$vcoefs; xcoefs0 = object$xcoefs; ucoefs = object$ucoefs; tcoefs = object$tcoefs
	tt = object$tms
	if (!inherits(object, "cgam")) { 
	        warning("calling predict.cgam(<fake-cgam-object>) ...")
        }
	if (missing(newData) | is.null(newData)) {
		etahat = object$etahat
		muhat = muhat.fun(etahat, fml = family$family)
		ans = list(fit = muhat, etahat = etahat, newbigmat = object$bigmat)
		return (ans) 
	}
	Terms = delete.response(tt)
	m = model.frame(Terms, newData)
	newdata = m
#new:
	newx0 = NULL; newxv = NULL; newx = NULL; newx_s = NULL; newu = NULL; newt = NULL; newz = NULL; newv = NULL
	rn = nrow(newdata)
	newetahat = NULL; newmuhat = NULL
	#newxbasis = matrix(nrow = nrow(newdata), ncol = (capm + capms)); 
	newxbasis = matrix(nrow = nrow(newdata), ncol = capm); 
	newubasis = matrix(nrow = nrow(newdata), ncol = capu); newtbasis = NULL; newbigmat = NULL
#######################
#local helper function#
#######################
	my_line = function(xp = NULL, y, x, end, start) {
		slope = NULL
		intercept = NULL
		yp = NULL
		slope = (y[end] - y[start]) / (x[end] - x[start])
		intercept = y[end] - slope * x[end]
		yp = intercept + slope * xp
		ans = new.env()
		ans$slope = slope
		ans$intercept = intercept
		ans$yp = yp
		ans
	}
	for (i in 1:ncol(newdata)) {
		if (is.null(attributes(newdata[,i])$shape)) {
#new:
			if (is.factor(newdata[,i])) {
				newdatai = as.numeric(levels(newdata[,i]))[newdata[,i]]
				if (length(levels(newdata[,i])) > 2) {
					klvls = length(levels(newdata[,i]))
					vals = as.numeric(levels(newdata[,i]))
					newimat = matrix(0, nrow = rn, ncol = klvls - 1)
					for (i1 in 1: (klvls - 1)) {
						for (i2 in 1:rn) {
							if (newdatai[i2] == vals[i1 + 1]) {
								newimat[i2, i1] = 1
							}
						}
					}
					newdatai = newimat
				}
				#newdatai = newimat
			} else {
				newdatai = newdata[,i]
			}
			newz = cbind(newz, newdatai)
			#newv = cbind(newv, newdatai)
		}
		if (is.numeric(attributes(newdata[,i])$shape)) {
			#if (attributes(newdata[,i])$shape != 17) {
				newx0 = cbind(newx0, newdata[,i])
			#}
			if ((attributes(newdata[,i])$shape > 2 & attributes(newdata[,i])$shape < 5) | (attributes(newdata[,i])$shape > 10 & attributes(newdata[,i])$shape < 13)) {
				#newv = cbind(newv, newdata[,i])
				newxv = cbind(newxv, newdata[,i])
			}
		} 
		if (is.character(attributes(newdata[,i])$shape)) {
       			if (attributes(newdata[,i])$shape == "tree") {
				newt = cbind(newt, newdata[,i])
			}
       			if (attributes(newdata[,i])$shape == "umbrella") {
				newu = cbind(newu, newdata[,i])	
			}
     		}
#new:
		#if (!is.null(newz)) {
		#	mat_cols <- ncol(newz)
		#	mat_rows <- nrow(newz)
		#	mat_rm <- NULL
		#	newz0 <- newz
		#	for (i in 1:mat_cols) {
		#		if (all(round(diff(newz0[, i]), 8) == 0)) {
		#			mat_rm <- c(mat_rm, i)
		 #     		}	
		  #  	}
		#	if (!is.null(mat_rm)) {
        	#		newz0 <- newz0[, -mat_rm, drop = FALSE]
		#	}
		#	newz <- newz0
		#}
		#newv <- cbind(newv, newz, newxv)

	}

#new:
if (!is.null(shapes)) {
	if (any(shapes == 17)) {
#print (shapes)
		kshapes <- length(shapes)
        	obs <- 1:kshapes
        	idx_s <- obs[which(shapes == 17)]; idx <- obs[which(shapes != 17)]
		
		newx1 <- newx0
		shapes0 <- 1:kshapes*0
 
		newx1[ ,1:length(idx_s)] <- newx0[ ,idx_s]
		shapes0[1:length(idx_s)] <- shapes[idx_s]
   
		if (length(idx) > 0) {
			newx1[ ,(1 + length(idx_s)):kshapes] <- newx0[ ,idx]
			shapes0[(1 + length(idx_s)):kshapes] <- shapes[idx]
    		}
		newx0 <- newx1; shapes <- shapes0
	}

#new code:

		if (all(shapes < 9)) {
			newx = newx0
			xid1 = xid10; xid2 = xid20
			xmat = xmat0
		} else if (all(shapes > 8)) {
			newx_s = newx0
			xid1_s = xid10; xid2_s = xid20
			xmat_s = xmat0
			numknots = numknots0
			knots = knots0
			sps = sps0
			sh = shapes
			ms = ms0
		} else if (any(shapes > 8) & any(shapes < 9)) {
			newx = newx0[, shapes < 9, drop = FALSE]; newx_s = newx0[, shapes > 8, drop = FALSE]
			xid1_s = xid10[shapes > 8]; xid2_s = xid20[shapes > 8]
			xid1 = xid10[shapes < 9]; xid2 = xid20[shapes < 9]
			xmat_s = xmat0[, shapes > 8, drop = FALSE]
			xmat = xmat0[, shapes < 9, drop = FALSE]
			numknots = numknots0[shapes > 8]
			knots = knots0[shapes > 8]
			sps = sps0[shapes > 8]
			sh = shapes[shapes > 8]
			ms = ms0[shapes > 8]
		}
}

#new:
	if (!is.null(newx_s)) {
 		ks = ncol(newx_s)
		delta = NULL; xs_coefs = NULL
		#xs_coefs = xcoefs0[shapes > 8] # check 
		for (i in 1:ks) {
			x0 = newx_s[ , i, drop = FALSE ]
			n0 = length(x0)
			pos1 = xid1_s[i] - (1 + capk + sum(shapes > 2 & shapes < 5 | shapes > 10 & shapes < 13))
			pos2 = xid2_s[i] - (1 + capk + sum(shapes > 2 & shapes < 5 | shapes > 10 & shapes < 13))
			xs_coefsi = xcoefs0[pos1:pos2]
			msi = ms[[i]]
			deltai_ans = makedelta(x0, sh[i], numknots[i], knots[[i]], space = sps[i], suppre = TRUE, interp = TRUE)
			deltai = matrix(deltai_ans$amat, ncol = n0)
			if (sh[i] > 10 & sh[i] < 13) {
				x = xmat_s[,i]
				xs = sort(x)
				ord = order(x)
				nx = length(x)
				obs = 1:nx
				m = nrow(msi)
				ms0 = matrix(0, nrow = m, ncol = n0)
				for (i1 in 1:n0) {
					for (i2 in 1:m) {
						ms0[i2, i1] = my_line(xp = x0[i1], y = msi[i2, ][ord], x = xs, end = nx, start = 1)$yp 
					}
				}
				deltai = deltai - ms0
			} else {
				deltai = deltai - msi
			}
			delta = rbind(delta, deltai)
			xs_coefs = c(xs_coefs, xs_coefsi)
		}
		etahat_s = t(delta) %*% xs_coefs
	}
#new:
	if (any(shapes == 17)) {
		vcoefs = vcoefs[1:(1 + capk + sum(shapes > 2 & shapes < 5 | shapes > 10 & shapes < 13))]
	}
#new:
	#if (!is.null(newv)) {
	#	mat_cols <- ncol(newv)
	#	mat_rows <- nrow(newv)
	#	mat_rm <- NULL
	#	newv0 <- newv
	#	for (i in 1:mat_cols) {
	#		if (all(round(diff(newv0[, i]), 8) == 0)) {
	#       	      mat_rm <- c(mat_rm, i)
	#      		}	
	#    	}
	#	if (!is.null(mat_rm)) {
        #		newv0 <- newv0[, -mat_rm, drop = FALSE]
	#	}
	#	newv <- newv0
	#}
	newv = cbind(1:rn*0 + 1, newz, newxv)
#print (newv)
#print (vcoefs)
	etahat = 1:rn*0
#######################################
#make newdata into same column as xmat#
#######################################
	#if (ncol(xmat) < 1) {
	#	stop ("There should be at least one constrained predictor to use predict.cgam!")
	#}
	if (!is.null(newx)) {
		newedge = NULL; xcoefs = NULL
		for (j in 1:ncol(xmat)) {
			x = xmat[,j]; nx = length(x); xs = sort(x)
			xu = unique(x); nxu = length(xu); xus = sort(xu)
			pos1 = xid1[j] - (1 + capk + sum(shapes > 2 & shapes < 5 | shapes > 10 & shapes < 13))
			pos2 = xid2[j] - (1 + capk + sum(shapes > 2 & shapes < 5 | shapes > 10 & shapes < 13))
			x_edges = t(bigmat[xid1[j]:xid2[j], , drop = FALSE])
			xcoefs = c(xcoefs, xcoefs0[pos1:pos2])
			#xcoefs = xcoefs0[pos1:pos2]
			nedge = ncol(x_edges)
			ord = order(x)
			ordu = order(xu)
			if (any(newx[,j] > max(x)) | any(newx[,j] < min(x))) {
				stop ("No extrapolation is allowed in cgam prediction!")
			}
			newedge0 = matrix(0, nrow = rn, ncol = nedge)
#######################################
#store the 2nd edge to check incr/decr#
#######################################
			if (nedge >= 2) {
				e2 = x_edges[, 2, drop = FALSE]; dist2 = round(diff(e2), 4); s2 = sign(dist2)
				e2ord = e2[ord]; dist2ord = round(diff(e2ord), 4); s2ord = sign(dist2ord)
			}
			for (i in 1:nedge) {	
				e = x_edges[, i, drop = FALSE]; dist = round(diff(e), 4); s = sign(dist)
				eord = e[ord]; distord = round(diff(eord), 4); sord = sign(distord)
				eu = unique(e); distu = round(diff(eu), 4); su = sign(distu)
				euord = eu[ordu]; distuord = round(diff(euord), 4); suord = sign(distuord)
				if (nx == 2) {
					newedge0[,i] = my_line(newx[,j], eord, xs, nx, 1)$yp
				} else {
###############################
#concave and convex: all edges#
###############################
					if (any(sord == 1) & any(sord == -1)) {
						for (l in 1:(nxu-2)) {
							if (suord[l] != suord[l+1]) {
								pos = l+1
								break
							}
						}
						if (any(newx[,j] >= xus[pos])) {
							ids = which(newx[,j] >= xus[pos])
							newedge0[ids,i] = my_line(newx[ids,j], y = euord, x = xus, end = nxu, start = pos)$yp
						} 
						if (any(newx[,j] < xus[pos])) {
							ids = which(newx[,j] < xus[pos])
							newedge0[ids,i] = my_line(newx[ids,j], y = euord, x = xus, end = pos, start = 1)$yp
						}  
#########################################################
#incr and decr:all edges except for the 1st and the last#
#########################################################
					} else if (sord[1] == 0 & sord[nx-1] == 0 & (sum(sord != 0) == 1)) {
						for (l in 1:(nx-2)) {
							if (sord[l] != sord[l+1]) {
								pos0 = l + 1
								pos = pos0 + 1
								break 
							} 
						}
						if (any(newx[,j] >= xs[pos])) {
							ids = which(newx[,j] >= xs[pos])
							newedge0[ids,i] = eord[pos]
						} 
						if (any(newx[,j] < xs[pos])) {
							ids = which(newx[,j] < xs[pos])
							newedge0[ids,i] = eord[pos0]
						}
#########################################################
#incr.conc,incr.conv,decr.conc, and decr.conv: all edges#
#########################################################
					} else { 
						pos0 = NULL
						for (l in 2:(nx-1)) {
							if (all(sord[l:(nx-1)] == 0) && any(sord[1:(l-1)] != 0)) {
								pos0 = l 
								break
							}
						}
						if (is.null(pos0)) {
							pos0 = NULL
							for (l in (nx-1):2) {
								if (all(sord[1:(l-1)] == 0) && any(sord[l:(nx-1)] != 0)) {
									pos0 = l
									break
								}
							}		
							if (is.null(pos0) & (all(suord > 0) | all(suord < 0))) {
								newedge0[,i] = my_line(newx[,j], y = eord, x = xs, end = nx, start = 1)$yp
							} else {
								obs1 = pos0:(nx-1)
								del = NULL
								if (any(sord[pos0:(nx-1)] == 0)) {
									id = obs1[sord[pos0:(nx-1)] == 0]
									del = c(del, id+1)	
								}
								if (is.null(del)) {
									pos = pos0
									if (any(newx[,j] >= xs[pos])) {
										ids = which(newx[,j] >= xs[pos])
										newedge0[ids,i] = my_line(newx[ids,j], y = eord, x = xs, end = nx, start = pos)$yp
									} 
									if (any(newx[,j] < xs[pos])) {
										ids = which(newx[,j] < xs[pos])
										newedge0[ids,i] = my_line(newx[ids,j], y = eord, x = xs, end = pos, start = 1)$yp
									} 
								} else {
									xc = unique(xs[-del])
									nxc = length(xc)
									ec = eord[-del]
									distc = round(diff(ec), 4)
									sc = sign(distc)
									for (l in 1:(nxc-2)) {
										if (sc[l] != sc[l+1]) {
											pos = l+1
											break
										}
									}
									if (any(newx[,j] >= xc[pos])) {
										ids = which(newx[,j] >= xc[pos])
										newedge0[ids,i] = my_line(newx[ids,j], y = ec, x = xc, end = nxc, start = pos)$yp
									} 
									if (any(newx[,j] < xc[pos])) {
										ids = which(newx[,j] < xc[pos])
										newedge0[ids,i] = my_line(newx[ids,j], y = ec, x = xc, end = pos, start = 1)$yp
									}
								}
							}
						} else {
							obs2 = 1:(pos0-1)
							del = NULL
							if (any(sord[1:(pos0-1)] == 0)) {
								id = obs2[sord[1:(pos0-1)] == 0]
								del = c(del, id+1)
							}
							if (is.null(del)) {
								pos = pos0
								if (any(newx[,j] >= xs[pos])) {
									ids = which(newx[,j] >= xs[pos])
									newedge0[ids,i] = my_line(newx[ids,j], y = eord, x = xs, end = nx, start = pos)$yp
								} 
								if (any(newx[,j] < xs[pos])) {
									ids = which(newx[,j] < xs[pos])
									newedge0[ids,i] = my_line(newx[ids,j], y = eord, x = xs, end = pos, start = 1)$yp
								} 
							} else {
								xc = unique(xs[-del])
								nxc = length(xc) 
								ec = eord[-del]
								distc = round(diff(ec), 4)
								sc = sign(distc)
								for (l in 1:(nxc-2)) {
									if (sc[l] != sc[l+1]) {
										pos = l+1
										break
									}
								}
								if (any(newx[,j] >= xc[pos])) {
									ids = which(newx[,j] >= xc[pos])
									newedge0[ids,i] = my_line(newx[ids,j], y = ec, x = xc, end = nxc, start = pos)$yp
								} 
								if (any(newx[,j] < xc[pos]))  {
									ids =  which(newx[,j] < xc[pos])
									newedge0[ids,i] = my_line(newx[ids,j], y = ec, x = xc, end = pos, start = 1)$yp
								}
							}
						}
					}
				}
			}		
##################################################################
#make corrections to the 1st and the last element of the newedge0# 
#for the 1st and the last edge of the incr and the decr case     #
##################################################################
			if (nx > 2 & nedge >= 2) {
				e1 = x_edges[, 1, drop = FALSE]; e1ord = e1[ord]
				ek = x_edges[, nedge, drop = FALSE]; ekord = ek[ord]
				pos01 = 1; pos1 = 2; pos0k = nx-1; posk = nx
				if (nx == 3 & nedge >= 2) { 
					if ((s2ord[1] != s2ord[nx-1]) & any(c(s2ord[1], s2ord[nx-1]) == 0)) {
						if (any(newx[,j] < xs[pos1] & newx[,j] >= xs[pos01])) {
							ids = which(newx[,j] < xs[pos1] & newx[,j] >= xs[pos01])
							newedge0[ids, 1] = e1ord[pos01]
						} 
						if (any(newx[,j] >= xs[pos0k] & newx[,j] < xs[posk])) {
							ids = which(newx[,j] >= xs[pos0k] & newx[,j] < xs[posk])
							newedge0[ids,nedge] = ekord[pos0k]
						}
					}
				} else if (nx > 3 & nedge >= 2) {
					if (s2ord[1] == 0 & s2ord[nx-1] == 0) {
						if (any(newx[,j] < xs[pos1] & newx[,j] >= xs[pos01])) {
							ids = which(newx[,j] < xs[pos1] & newx[,j] >= xs[pos01])
							newedge0[ids,1] = e1ord[pos01]
						} 
				 		if (any(newx[,j] >= xs[pos0k] & newx[,j] < xs[posk])) {
							ids = which(newx[,j] >= xs[pos0k] & newx[,j] < xs[posk])
							newedge0[ids,nedge] = ekord[pos0k]
						}
					}
				}
			}	
		newedge = cbind(newedge, newedge0)
		#newxbasis[, (xid1[j] - (np - capms)):(xid2[j] - (np - capms))] = newedge0
		#etahat = etahat + newedge0 %*% xcoefs 
		}
#print (dim(cbind(newv, newedge)))
#print (length(c(vcoefs, xcoefs)))
#print (sum(is.na(etahat)))
		#etahat = newv %*% vcoefs + etahat
		etahat = cbind(newv, newedge) %*% c(vcoefs, xcoefs) 
	}
	if (!is.null(newu)) {
		newuedge = NULL
		for (j in 1:ncol(umb)) {
			u = umb[,j]
			nu = length(u)
			us = sort(u)
			ord = order(u)
			if (any(newu[,j] > max(u)) | any(newu[,j] < min(u))) {
				stop ("no extrapolation is allowed in cgam!")
			}
			#u_edges = t(umbrella.fun(u))
			pos1 = uid1[j]; pos2 = uid2[j]
			u_edges = t(bigmat[pos1:pos2, , drop = FALSE])
			nuedge = ncol(u_edges)
			newuedge0 = matrix(0, nrow = rn, ncol = nuedge)
			for (i in 1:nuedge) {
				ue = u_edges[, i, drop = FALSE]; udist = round(diff(ue), 4); s = sign(udist)
				ueord = ue[ord]; udistord = round(diff(ueord), 4); sord = sign(udistord)
				obs = 1:(nu-1)
				posin = obs[sord != 0] + 1
				pos = unique(c(1, posin, nu))
				uk = us[pos]
				uek = ueord[pos]
				npos = length(pos)
				nd = length(newu[,j])
				for (l in 1:(npos-1)) {
					if (any(newu[,j] == uk[npos])) {
						ids = which(newu[,j] == uk[npos])
						newuedge0[ids,i] = uek[npos]
					}
					if (any(newu[,j] >= uk[l]) & any(newu[,j] < uk[l+1])) {
						ids = which(newu[,j] >= uk[l] & newu[,j] < uk[l+1])
						newuedge0[ids,i] = uek[l]
					}
				}
			}
		newuedge = cbind(newuedge, newuedge0)
		newubasis[, (uid1[j] - np - capm):(uid2[j] - np - capm)] = newuedge0 
		}
		etahat = etahat + newuedge %*% ucoefs
	}
	if (!is.null(newt)) {
		tru = unique(tr); placebo = min(tru)
		pos_end = 0
		tr_etahat = 0
		for (i in 1: ncol(newt)) {
			#te = t(tree.fun(tr[,i]))
			pos1 = tid1[i]; pos2 = tid2[i]
			te = t(bigmat[pos1:pos2, , drop = FALSE])
			pos_start = pos_end + 1
			pos_end = pos_end + ncol(te)
			newtu = unique(newt[,i])
			for (inewt in 1: length(newtu)) {
				if (!any(tru[,i] == newtu[inewt])) {
					stop ("new tree ordering factor must be among the old tree ordering factors!")
				} 
			} 
			if (any(newt[,i] != placebo)) {	
				tr_k = sum(newtu != placebo)				
				tr_etahat = tr_etahat + te[, 1:tr_k , drop = FALSE] %*% c(tcoefs[pos_start:pos_end][1:tr_k])
			}
		}
		etahat = etahat + tr_etahat
	} 
	#newetahat = c(newetahat, etahat)
	newetahat = etahat 	
	if (is.null(newx)) {
		newetahat = newv %*% vcoefs + newetahat
	}
	if (!is.null(newx_s)) {
		newetahat = etahat_s + newetahat
	}
	newmuhat = muhat.fun(newetahat, fml = family$family)
	#if (!is.null(newt)) {
	#	newtbasis = t(tree.delta[ ,1:nrow(newData), drop = FALSE])
	#}
	#ans = new.env()
	#if (!is.null(newt)) {
	#	newbigmat = t(cbind(newv, newxbasis, newubasis, newtbasis))
	#} else {
	#	newbigmat = t(cbind(newv, newxbasis, newubasis))
	#}
	#ans = list(newv = newv, newxbasis = newxbasis, newubasis = newubasis, newtbasis = newtbasis, newbigmat = newbigmat, etahat = newetahat, muhat = newmuhat)
	ans = list(etahat = newetahat, muhat = newmuhat)
	return (ans) 
}

###########################################
#create a 3D plot for a cgam or wps object#
###########################################
plotpersp <- function(object, x1, x2,...)UseMethod("plotpersp")

################
#plotpersp.cgam#
################ 
plotpersp.cgam <- function(object, x1 = NULL, x2 = NULL, data = NULL, surface = "mu", categ = NULL, col = NULL, random = FALSE, x_grid = 10, y_grid = 10, xlim = range(x1), ylim = range(x2), zlim = NULL, xlab = NULL, ylab = NULL, zlab = NULL, main = NULL, sub = NULL, th = NULL, phi = 15, r = sqrt(3), d = 1, scale = TRUE, expand = 1, border = NULL, ltheta = NULL, lphi = 0, shade = NA, box = TRUE, axes = TRUE, nticks = 5, ticktype = "detailed",...) {
	if (!inherits(object, "cgam")) { 
	        warning("calling plotpersp(<fake-cgam-object>) ...")
        }
#new:
#if (is.null(x1) | is.null(x2)) {
#	stop ("x1 and x2 must be provided!")
#}
	cl <- match.call()
	nms <- cl[-c(1, 2)]
	#lnms <- length(nms)
	x1nm <- nms[1]$x
	x1nm <- deparse(x1nm)
	x2nm <- nms[2]$x
	x2nm <- deparse(x2nm)
	ynm <- object$ynm
	xmat <- object$xmat
	is_fac <- object$is_fac
	is_param <- object$is_param
	family <- object$family
	fml <- family$family
	cicfamily <- CicFamily(family)
	muhat.fun <- cicfamily$muhat.fun
	znms <- object$znms
	kznms <- length(znms)
#print (class(categ))
	if (!is.null(categ)) {
		if (!is.character(categ)) {
			warning("categ must be a character argument!")
		} else if (!any(znms == categ)) {
			warning(paste(categ, "is not an exact character name defined in the cgam fit!"))
		}
	}
#new: no layer for linear term
if (!is.null(categ)) {
	obsz = 1:kznms
	zid = obsz[znms == categ]
	if (!(is_fac[zid])) {
		categ = NULL
	}
}
	shapes <- object$shapes
#new:	
	#zid1 = object$zid1 - 1 - length(shapes)
	#zid2 = object$zid2 - 1 - length(shapes)
	zid1 <- object$zid1
	zid2 <- object$zid2

	kznms <- length(znms)
	zmat <- object$zmat
	#zcoefs = object$zcoefs
#new: exclude the coef for the one vector
	zcoefs <- object$zcoefs[-1]
	xmat0 <- object$xmat
	xmatnms <- object$xmatnms
	knms <- length(xmatnms)	
	obs <- 1:knms

	#if (!any(xmatnms == x1nm)) {
	#	warning(paste(x1nm, "is not an exact character name defined in the cgam fit!"))
	#}
	#if (!any(xmatnms == x2nm)) {
	#	warning(paste(x2nm, "is not an exact character name defined in the cgam fit!"))
	#}
	#x1id = obs[xmatnms == x1nm]
	#x2id = obs[xmatnms == x2nm]
	if (!is.null(data)) {
		if (!is.data.frame(data)) {
			stop ("User need to make the data argument a data frame with names for each variable!")
		}
		datnms <- names(data)
		if (!any(datnms == x1nm) | !any(datnms == x2nm)) {
			stop ("Check the accuracy of the names of x1 and x2!")
		}
		x1 <- data[ ,which(datnms == x1nm)]
		x2 <- data[ ,which(datnms == x2nm)]
		bool <- apply(xmat, 2, function(x) all(x1 == x))
		if (any(bool)) {
			x1id <- obs[bool]
		}
		bool <- apply(xmat, 2, function(x) all(x2 == x))
		if (any(bool)) {
			x2id <- obs[bool]
		}
	} else {
		#if (any(xmatnms == x1nm)) {
		#	x1id <- obs[xmatnms == x1nm]
		#} else {
		#	bool <- apply(xmat, 2, function(x) all(x1 == x))
		#	if (any(bool)) {
		#		x1id <- obs[bool]
		#	}
		#}
		#if (any(xmatnms == x2nm)) {
		#	x2id <- obs[xmatnms == x2nm]
		#} else {
		#	bool <- apply(xmat, 2, function(x) all(x2 == x))
		#	if (any(bool)) {
		#		x2id <- obs[bool]
		#	}
		#}
#new: x1 and x2 are in .Globe, not in formula
		if (all(xmatnms != x1nm)) {
			bool <- apply(xmat, 2, function(x) all(x1 == x))
			if (any(bool)) {
				x1id <- obs[bool]
#change x1nm to be the one in formula
				x1nm <- xmatnms[bool] 
			} else {
				stop (paste(paste("'", x1nm, "'", sep = ''), "is not a predictor defined in the cgam fit!"))
			}
		} else {
			x1id <- obs[xmatnms == x1nm]
		}
		if (all(xmatnms != x2nm)) {
			bool <- apply(xmat, 2, function(x) all(x2 == x))
			if (any(bool)) {
				x2id <- obs[bool]
				x2nm <- xmatnms[bool] 
			} else {
				stop (paste(paste("'", x2nm, "'", sep = ''), "is not a predictor defined in the cgam fit!"))
			}
		} else {
			x2id <- obs[xmatnms == x2nm]
		}
	}
	xmat <- cbind(x1, x2)
	x1g <- 0:x_grid / x_grid * .95 * (max(xmat[,1]) - min(xmat[,1])) + min(xmat[,1]) + .025 * (max(xmat[,1]) - min(xmat[,1]))
 	n1 <- length(x1g)
	x2g <- 0:y_grid / y_grid * .95 * (max(xmat[,2]) - min(xmat[,2])) + min(xmat[,2]) + .025 * (max(xmat[,2]) - min(xmat[,2]))
	n2 <- length(x2g)	

	xgmat <- matrix(nrow = n1, ncol = n2)
	eta0 <- object$coefs[1]
	thvecs <- object$etacomps
	for (i2 in 1:n2) {
		for (i1 in 1:n1) {
			x1a <- max(xmat[xmat[,1] <= x1g[i1], 1])
			x1b <- min(xmat[xmat[,1] > x1g[i1], 1])
			v1a <- min(thvecs[x1id, xmat[,1] == x1a])
			v1b <- min(thvecs[x1id, xmat[,1] == x1b])
			alp <- (x1g[i1] - x1a) / (x1b - x1a)
			th1add <- (1 - alp) * v1a + alp * v1b
			x2a <- max(xmat[xmat[,2] <= x2g[i2],2])
			x2b <- min(xmat[xmat[,2] > x2g[i2],2])
			v2a <- min(thvecs[x2id, xmat[,2] == x2a])
			v2b <- min(thvecs[x2id, xmat[,2] == x2b])
			alp <- (x2g[i2] - x2a) / (x2b - x2a)
			th2add <- (1 - alp) * v2a + alp * v2b	
			xgmat[i1,i2] <- eta0 + th1add + th2add
		}
	}
	x3_add <- 0
	if (knms >= 3) {
		x3id <- obs[-c(x1id, x2id)]
		kx3 <- length(x3id)
		for (i in 1:kx3) {
			x3i <- xmat0[, x3id[i]]
			x3i_use <- max(x3i[x3i <= median(x3i)])
			x3i_add <- min(thvecs[x3id[i], x3i == x3i_use])			
			x3_add <- x3_add + x3i_add
		}
	} 
	if (surface == "eta") {
		xgmat <- xgmat + as.numeric(x3_add)
	}	
	if (is.null(categ) & surface == "mu") {
		z_add <- 0
		if (!is.null(znms)) {
			kzids <- length(zid1)
			for (i in 1:kzids) {
				pos1 <- zid1[i]; pos2 <- zid2[i]
				zi <- zmat[, pos1:pos2, drop = FALSE]
				zcoefsi <- zcoefs[pos1:pos2]
				for (j in 1:ncol(zi)){
					uzij <- unique(zi[,j])
					kuzij <- length(uzij)
					nmodej <- sum(zi[,j] == uzij[1])
					zij_mode <- uzij[1]
					for (u in 2:kuzij) {
						if (sum(zi[,j] == uzij[u]) > nmodej) {
							zij_mode <- uzij[u]
							nmodej <- sum(zi[,j] == uzij[u])
						}
					}
					obsuzij <- 1:length(uzij)
					uzhatij <- uzij * zcoefsi[j] 
					zij_add <- uzhatij[obsuzij[uzij == zij_mode]]	
					z_add <- z_add + zij_add 
				}
			}
		}
		xgmat <- xgmat + as.numeric(x3_add) + as.numeric(z_add)
		for (i2 in 1:n2) {
			for (i1 in 1:n1) {
				xgmat[i1, i2] <- muhat.fun(xgmat[i1, i2], fml = fml)	
			}
		}
	} else if (!is.null(categ) & surface == "mu"){
		xgmats <- list()
		mins <- NULL; maxs <- NULL
		obsz <- 1:kznms
		zid <- obsz[znms == categ]
#print (class(znms[znms == categ]))
		pos1 <- zid1[zid]; pos2 <- zid2[zid]
#print (pos1)
#print (pos2)
		zi <- zmat[, pos1:pos2, drop = FALSE]
		z_add <- 1:nrow(zi)*0
		zcoefsi <- zcoefs[pos1:pos2]
#print (zcoefsi)
		for (j in 1:ncol(zi)) {
			zij <- zi[,j]
			zijhat <- zij * zcoefsi[j]
			z_add <- z_add + zijhat
		}
		z_add <- unique(z_add)
		kz_add <- length(z_add)
#new: plot the smallest one first:
		z_add <- z_add[order(z_add)]
#print (z_add)
		for (iz in 1:kz_add) {
			xgmats[[iz]] <- xgmat + as.numeric(x3_add) + z_add[iz]
			mins <- c(mins, min(xgmats[[iz]]))
			maxs <- c(maxs, max(xgmats[[iz]]))
			for (i2 in 1:n2) {
				for (i1 in 1:n1) {
					xgmats[[iz]][i1, i2] <- muhat.fun(xgmats[[iz]][i1, i2], fml = fml)	
				}
			}

		}
	}
	if (is.null(xlab)) {
		#xlab = deparse(x1nm)
		xlab <- x1nm
	}
	if (is.null(ylab)) {
		#ylab = deparse(x2nm)
		ylab <- x2nm
	}
	if (is.null(zlab)) {
		if (surface == "mu") {
			if (fml == "binomial") {
				zlab <- paste("Pr(", ynm, ")")
			} else if (fml == "poisson" | fml == "gaussian") {
				zlab <- paste("est mean of", ynm)
			}		
		}
		if (surface == "eta") {
			if (fml == "binomial") {
				zlab <- paste("est log odds ratio of", ynm)
			}  else if (fml == "poisson") {
				zlab <- paste("est log mean of", ynm)
			} else if (fml == "gaussian") {
				zlab <- paste("est mean of", ynm)
			}
		}
	}
	if (is.null(zlim)) {
		zlim <- range(xgmat, na.rm = TRUE)
	}
	palette <- c("peachpuff", "lightblue", "limegreen", "grey", "wheat", "yellowgreen", "seagreen1", "palegreen", "azure", "whitesmoke")
	if (!is.null(categ) & surface == "mu") {
		#palette = c("peachpuff", "lightblue", "grey", "wheat", "yellowgreen", "plum", "limegreen", "paleturqoise", "azure", "whitesmoke")
		kxgm <- length(xgmats)
		if (is.null(col)) {
			#if (kxgm == 2) {
			#	col = c("peachpuff", "lightblue")
			#} else if (kxgm == 3) { 
			#	col = c("peachpuff", "lightblue", "grey")
			#} else if (kxgm > 3 & kxgm < 11) {
			#	col = sample(palette, replace = FALSE)
			if (random) {
				col <- sample(palette, size = kxgm, replace = FALSE)
#print (col)
			} else {
#print (kxgm)
				if (kxgm > 1 & kxgm < 11) {
					col <- palette[1:kxgm]
				} else {
					integ <- floor(kxgm / 10)
					rem <- kxgm %% 10
					kint <- length(integ)
					#col = character(length = kxgm)
#print (col)
					col <- NULL
					for (i in 1:kint) {
#print (col[1 + (i - 1) * 10: i * 10])
#print (palette)
						#col[1 + (i - 1) * 10: i * 10] = palette
						col <- c(col, palette)
					}
#print (col)
					#col[(kint * 10 + 1):kxgm] = palette[(kint * 10 + 1):kxgm]
					col <- c(col, palette[1:rem])
#print ((kint * 10 + 1):kxgm)
#print (col)
#print (integ)
				}
			} 
		} else {
			if (length(col) < kxgm) {
				rem <- kxgm - length(col)
				nrem <- length(rem)
				rem_col <- palette[1:nrem]
				col <- c(col, rem_col)
			} else if (length(col) > kxgm) {
				col <- col[1:kxgm]
				print (paste("The first", kxgm, "colors are used!"))
			}
			if (random) {
				print ("User defined colors are used!")
			}
		}
#print (kxgm)
		for (i in 1:kxgm) {
			xgmat <- xgmats[[i]]
			if (is.null(th) | !is.numeric(th)) {
				th <- -40
			}
			if (is.null(ltheta) | !is.numeric(ltheta)) {
				ltheta <- -135
			}
			persp(x1g, x2g, xgmat, xlim = xlim, ylim = ylim, zlim = c(min(mins), max(maxs)), xlab = xlab, ylab = ylab, zlab = zlab, main = main, sub = sub, theta = th, phi = phi, r = r, d = d, scale = scale, expand = expand, col = col[i], border = border, ltheta = ltheta, lphi = lphi, shade = shade, box = box, axes = axes, nticks = nticks, ticktype = ticktype)
			par(new = TRUE)
		}
	par(new = FALSE)
	} else {
		if (is.null(col)) {
			if (random) {
				col <- sample(palette, size = 1, replace = FALSE)
			} else {
				col <- "white"
			}
		} else {
			if (length(col) > 1) {
				col <- col[1]
				print ("The first color is used!")
			} 
			if (random) {
				print ("User defined color is used!")
			}
		}
		if (is.null(th) | !is.numeric(th)) {
			th <- -40
		}
		if (is.null(ltheta) | !is.numeric(ltheta)) {
			ltheta <- -135
		}
		persp(x1g, x2g, xgmat, xlim = xlim, ylim = ylim, zlim = zlim, xlab = xlab, ylab = ylab, zlab = zlab, main = main, sub = sub, theta = th, phi = phi, r = r, d = d, scale = scale, expand = expand, col = col, border = border, ltheta = ltheta, lphi = lphi, shade = shade, box = box, axes = axes, nticks = nticks, ticktype = ticktype)
	}
#print (col)
}

#################
#plotpersp.wps#
################ 
plotpersp.wps = function(object, x1 = NULL, x2 = NULL, data = NULL, surface = "C", categ = NULL, col = NULL, random = FALSE, xlim = range(x1), ylim = range(x2), zlim = NULL, xlab = NULL, ylab = NULL, zlab = NULL, main = NULL, sub = NULL, th = NULL, phi = 15, r = sqrt(3), d = 1, scale = TRUE, expand = 1, border = NULL, ltheta = NULL, lphi = 0, shade = NA, box = TRUE, axes = TRUE, nticks = 5, ticktype = "detailed",...) {
	if (!inherits(object, "wps")) { 
	        warning("calling plotpersp(<fake-wps-object>) ...")
        }
#new:
#if (is.null(x1) | is.null(x2)) {
#	stop ("x1 and x2 must be provided!")
#}
	cl = match.call()
	nms = cl[-c(1, 2)]
	#lnms = length(nms)
	xnms = object$xnms
	#knms = length(xnms)	
	#obs = 1:knms
	x1nm0 = nms[1]$x
	x1nm0 = deparse(x1nm0)
	x2nm0 = nms[2]$x
	x2nm0 = deparse(x2nm0)
	xmat = object$xmat
	is_fac = object$is_fac
	#is_param = object$is_param
#new: switch xnms
if (!is.null(data)) {
	if (!is.data.frame(data)) {
		stop ("User need to make the data argument a data frame with names for each variable!")
	}
	datnms = names(data)
	if (!any(datnms == x1nm0) | !any(datnms == x2nm0)) {
		stop ("Check the accuracy of the names of x1 and x2!")
	}
	x1 = data[ ,which(datnms == x1nm0)]
	x2 = data[ ,which(datnms == x2nm0)]
	#if (x1nm0 != xnms[1] & x2nm0 != xnms[2]) {
	#	x1nm = x2nm0
	#	x2nm = x1nm0
	#	tmp = x1
	#	x1 = x2
	#	x2 = tmp
	#} else {x1nm = x1nm0; x2nm = x2nm0}
} else {
	if (all(xnms != x1nm0)) {
		#stop (paste(paste("'", x1nm0, "'", sep = ''), "is not an exact predictor name defined in the cgam fit!"))
		bool = apply(xmat, 2, function(x) all(x1 == x))
		if (any(bool)) {
			#x1id = obs[bool]
			x1nm0 = xnms[bool] 
		} else {
			stop (paste(paste("'", x1nm0, "'", sep = ''), "is not a predictor defined in the wps fit!"))
		}
	} 
	if (all(xnms != x2nm0)) {
		#stop (paste(paste("'", x2nm0, "'", sep = ''), "is not an exact predictor name defined in the cgam fit!"))
		bool = apply(xmat, 2, function(x) all(x2 == x))
		if (any(bool)) {
			#x2id = obs[bool]
			x2nm0 = xnms[bool] 
		} else {
			stop (paste(paste("'", x2nm0, "'", sep = ''), "is not a predictor defined in the wps fit!"))
		}
	} 
} 
if (x1nm0 != xnms[1] &  x2nm0 != xnms[2]) {
	x1nm = x2nm0
	x2nm = x1nm0
	tmp = x1
	x1 = x2
	x2 = tmp
} else {x1nm = x1nm0; x2nm = x2nm0}
#print (x1nm0)
#print (x2nm0)
	ynm = object$ynm
#xmat is delta
	#delta = object$delta
	znms = object$znms
	decrs = object$decrs
	kznms = length(znms)
#zmat include 1 vector
	zmat = object$zmat
	zmat0 = zmat[, -1, drop = FALSE]
	zcoefs = object$zcoefs[-1]
	zid1 = object$zid1
	zid2 = object$zid2
	ah = object$coefs
	ahu = object$coefsu
	k1 = object$k1
	k2 = object$k2
	m1 = length(k1)
	m2 = length(k2)
	p = dim(zmat)[2]
	apl = 1:(m1 + m2 - 1 + (m1 - 1) * (m2 - 1))
	aplu = 1:(m1 + m2 - 1 + (m1 - 1) * (m2 - 1))
	apl[1] = ah[1]
	aplu[1] = ahu[1]
	apl[2:(m1 + m2 - 1 + (m1 - 1) * (m2 - 1))] = ah[(p + 1):(m1 + m2 - 1 + (m1 - 1) * (m2 - 1) + p - 1)]
	aplu[2:(m1 + m2 - 1 + (m1 - 1) * (m2 - 1))] = ahu[(p + 1):(m1 + m2 - 1 + (m1 - 1) * (m2 - 1) + p - 1)]
	mupl = matrix(apl[1], nrow = m1, ncol = m2)
	muplu = matrix(aplu[1], nrow = m1, ncol = m2)
	for (i1 in 2:m1) {
		mupl[i1, ] = mupl[i1, ] + apl[i1]
		muplu[i1, ] = muplu[i1, ] + aplu[i1]
	}
	for (i2 in 2:m2) {
		mupl[, i2] = mupl[, i2] + apl[m1 - 1 + i2]
		muplu[, i2] = muplu[, i2] + aplu[m1 - 1 + i2]
	}
	for (i1 in 2:m1) {
		for (i2 in 2:m2) {
			mupl[i1, i2] = mupl[i1, i2] + apl[m1 + m2 - 2 + (i1 - 2) * (m2 - 1) + i2]
			muplu[i1, i2] = muplu[i1, i2] + aplu[m1 + m2 - 2 + (i1 - 2) * (m2 - 1) + i2]
		}
	}
## reverse transform for decreasing
	if (is.null(th)) {
		ang = NULL
		if (!decrs[1] & !decrs[2]) {
			if (is.null(ang)) {
				ang = -40
			}
		} else if (decrs[1] & !decrs[2]) { 
			#k1 = -k1[m1:1]; 
			mupl = mupl[m1:1, 1:m2]; muplu = muplu[m1:1, 1:m2]
			if (is.null(ang)) {
				ang = 40
			}
		} else if (!decrs[1] & decrs[2]) {
			#k2 = -k2[m2:1]; 
			mupl = mupl[1:m1, m2:1]; muplu = muplu[1:m1, m2:1]
			if (is.null(ang)) {
				ang = 240
			}
		} else if (decrs[1] & decrs[2]) {
			#k1 = -k1[m1:1]; k2 = -k2[m2:1]; 
			mupl = mupl[m1:1, m2:1]; muplu = muplu[m1:1, m2:1]
			if (is.null(ang)) {
				ang = 140
			}
		}
	} else {ang = th}
#zmat: 
#new: no layer for linear term
if (!is.null(categ)) {
	obsz = 1:kznms
	zid = obsz[znms == categ]
	if (!(is_fac[zid])) {
		categ = NULL
	}
}
	if (is.null(categ)) {
		z_add = 0
		if (!is.null(znms)) {
			kzids = length(zid1)
			for (i in 1:kzids) {
				pos1 = zid1[i]; pos2 = zid2[i]
#zi is a factor 
				zi = zmat0[, pos1:pos2, drop = FALSE]
				zcoefsi = zcoefs[pos1:pos2]
				for (j in 1:ncol(zi)){
#find the 'mode' of the jth column of zi; add the coef corresponding to the 'mode'
					uzij = unique(zi[,j])
					kuzij = length(uzij)
					nmodej = sum(zi[,j] == uzij[1])
					zij_mode = uzij[1]
					for (u in 2:kuzij) {
						if (sum(zi[,j] == uzij[u]) > nmodej) {
							zij_mode = uzij[u]
							nmodej = sum(zi[,j] == uzij[u])
						}
					}
					obsuzij = 1:length(uzij)
					uzhatij = uzij * zcoefsi[j] 
					zij_add = uzhatij[obsuzij[uzij == zij_mode]]	
					z_add = z_add + zij_add 
				}
			}
		}
		mupl = mupl + as.numeric(z_add)
		muplu = muplu + as.numeric(z_add)
		mins = min(mupl); maxs = max(mupl)
		minsu = min(muplu); maxsu = max(muplu)
	} else {
		mupls = muplus = list()
		mins = maxs = NULL
		minsu = maxsu = NULL		
		obsz = 1:kznms
		zid = obsz[znms == categ]
		pos1 = zid1[zid]; pos2 = zid2[zid]
		zcoefsi = zcoefs[pos1:pos2]
#include the base level
		zcoefsi = c(0, zcoefsi)
		z_add = sort(zcoefsi)
		kz_add = length(z_add)
		for (iz in 1:kz_add) {
			mupls[[iz]] = mupl + z_add[iz]
			mins = c(mins, min(mupls[[iz]]))
			maxs = c(maxs, max(mupls[[iz]]))		
			muplus[[iz]] = muplu + z_add[iz]
			minsu = c(minsu, min(muplus[[iz]]))
			maxsu = c(maxsu, max(muplus[[iz]]))
		}

	}
	palette = c("peachpuff", "lightblue", "limegreen", "grey", "wheat", "yellowgreen", "seagreen1", "palegreen", "azure", "whitesmoke")
	if (is.null(categ)) {
		if (is.null(col)) {
			if (random) {
				col = sample(palette, size = 1, replace = FALSE)
			} else {
				col = "white"
			}
		} 
		if (surface == 'C') {
			musurf = mupl
			if (is.null(main)) {
				main = 'Constrained Warped-Plane Spline Surface'
			}
			zlim0 = c(min(mins), max(maxs))
		} else if (surface == 'U') {
			musurf = muplu
			if (is.null(main)) {
				main =  'Unconstrained Warped-Plane Spline Surface'
			}
			zlim0 = c(min(minsu), max(maxsu))
		}
		persp(k1, k2, musurf, xlim = xlim, ylim = ylim, zlim = zlim0, xlab = x1nm, ylab = x2nm, zlab = ynm, main = main,
sub = sub, theta = ang, phi = phi, r = r, d = d, scale = scale, col = col, border = border, ltheta = ltheta, lphi = lphi, shade = shade, box = box, axes = axes, nticks = nticks, ticktype = ticktype, cex.axis = .75, font.main = 1, expand = 1)
	} else {
		kxgm = length(mupls)
		if (is.null(col)) {
			if (random) {
				col = sample(palette, size = kxgm, replace = FALSE)
			} else {
				if (kxgm > 1 & kxgm < 11) {
					col = palette[1:kxgm]
				} else {
					integ = floor(kxgm / 10)
					rem = kxgm %% 10
					kint = length(integ)
					col = NULL
					for (i in 1:kint) {
						col = c(col, palette)
					}
					col = c(col, palette[1:rem])
				}
			} 
		} else {
			if (length(col) < kxgm) {
				rem = kxgm - length(col)
				nrem = length(rem)
				rem_col = palette[1:nrem]
				col = c(col, rem_col)
			} else if (length(col) > kxgm) {
				col = col[1:kxgm]
				print (paste("The first", kxgm, "colors are used!"))
			}
			if (random) {
				print ("User defined colors are used!")
			}
		}
		for (i in 1:kxgm) {
			mupli = mupls[[i]]
			muplui = muplus[[i]]
			if (surface == 'C') {
				musurf = mupli
				if (is.null(main)) {
					main = 'Constrained Warped-Plane Spline Surface'
				}
				zlim0 = c(min(mins), max(maxs))
			} else if (surface == 'U') {
				musurf = muplui
				if (is.null(main)) {
					main = 'Unconstrained Warped-Plane Spline Surface' 
				}				
				zlim0 = c(min(minsu), max(maxsu))
			}
			par(mar = c(4, 2, 2, 2))
			persp(k1, k2, musurf, xlim = xlim, ylim = ylim, zlim = zlim0, xlab = x1nm, ylab = x2nm, zlab = ynm, main = main,
sub = sub, theta = ang, phi = phi, r = r, d = d, scale = scale, col = col[i], border = border, ltheta = ltheta, lphi = lphi, shade = shade, box = box, axes = axes, nticks = nticks, ticktype = ticktype, cex.axis = .75, font.main = 1, expand = 1)
			par(new = TRUE)
		}
		par(new = FALSE)
	}
}


#####
#wps#
#####
wps <- function(formula, family = gaussian(), data = NULL, weights = NULL, pnt = FALSE, pen = 0, cpar = 1.5)
{
  cl <- match.call()
  if (is.character(family)) 
     family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) 
     family <- family()
  if (is.null(family$family)) 
     stop("'family' not recognized!")
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  ynm <- names(mf)[1]
  mt <- attr(mf, "terms")
  y <- model.response(mf, "any")
  if (family$family == "binomial") {
	if (class(y) == "factor") {
		y = ifelse(y == levels(y)[1], 0, 1)
	}
  }
  warp.delta <- NULL
  nums <- NULL; ks <- list(); sps <- NULL
  xmat <- NULL; x1 <- NULL; x2 <- NULL; xnms <- NULL; nks0 <- NULL; ks0 <- NULL; dcs <- NULL
  zmat <- NULL; zid <- NULL; zid1 <- NULL; zid2 <- NULL; znms <- NULL; is_fac <- NULL; is_param <- NULL; vals <- NULL; st <- 1; ed <- 1
  for (i in 2: ncol(mf)) {
    if (is.character(attributes(mf[, i])$shape)) {
	sps <- attributes(mf[, i])$space
	dcs <- attributes(mf[, i])$decreasing
	ks0 <- attributes(mf[, i])$knots
	nks0 <- attributes(mf[, i])$numknots
	x1 <- (mf[, i])[, 1]
	x2 <- (mf[, i])[, 2]
	xmat <- cbind(x1, x2)
	xnms <- c(xnms, attributes(mf[, i])$name)
	ans_warp <- makedelta_wps(x1t = x1, x2t = x2, m1_0 = nks0[1], m2_0 = nks0[2], k1 = ks0$k1, k2 = ks0$k2, space = sps, decreasing = dcs)
	warp.delta0 <- ans_warp$delta
	k1 <- ans_warp$k1
	k2 <- ans_warp$k2
	ks[[1]] <- k1
	ks[[2]] <- k2
	warp.delta <- cbind(warp.delta, warp.delta0)
    }
    if (is.null(attributes(mf[, i])$shape)) {
	if (!is.null(names(mf)[i])) {
	  znms <- c(znms, names(mf)[i])
	}
        if (!is.matrix(mf[, i])) {
      		zid <- c(zid, i)
		is_param <- c(is_param, TRUE)
        	if (is.factor(mf[, i])) {
	  		is_fac <- c(is_fac, TRUE)
	 	 	ch_char <- suppressWarnings(is.na(as.numeric(levels(mf[, i]))))
          		if (any(ch_char)) {
	    			vals <- c(vals, unique(levels(mf[, i]))[-1])
          		} else {
	    			vals <- c(vals, as.numeric(levels(mf[, i]))[-1])
	  		}
          		nlvs <- length(attributes(mf[, i])$levels)
			ed <- st + nlvs - 2 
			zid1 <- c(zid1, st)
	  		zid2 <- c(zid2, ed)
	  		st <- st + nlvs - 1
	  		zmat0 <- model.matrix(~ mf[, i])[, -1, drop = FALSE]
	  		zmat <- cbind(zmat, zmat0)
    		} else {
			is_fac <- c(is_fac, FALSE)
			zmat <- cbind(zmat, mf[, i])
			ed <- st
			zid1 <- c(zid1, st)
			zid2 <- c(zid2, ed)	
			st <- st + 1		
			vals <- c(vals, "")
            	}
	} else {
	 	is_param <- c(is_param, FALSE)
		is_fac <- c(is_fac, FALSE)
	  	zmat0 <- mf[, i]
	  	zmat <- cbind(zmat, zmat0)
	  	#vals <- c(vals, 1)
	  	zid <- c(zid, i)
	  	nlvs <- ncol(zmat0) + 1
	  	ed <- st + nlvs - 2
	  	zid1 <- c(zid1, st)
	  	zid2 <- c(zid2, ed)
	  	st <- st + nlvs - 1
	}
    }
  }
  ans <- wps.fit(x1t = x1, x2t = x2, y = y, zmat = zmat, w = weights, pen = pen, pnt = pnt, cpar = cpar, decrs = dcs, delta = warp.delta, kts = ks)
  rslt <- list(k1 = ans$k1, k2 = ans$k2, muhat = ans$muhat, muhatu = ans$muhatu, SSE1 = ans$sse1, SSE0 = ans$sse0, edf = ans$edf, edfu = ans$edfu, delta = ans$delta, zmat = ans$zmat, xmat = xmat, coefs = ans$coefs, coefsu = ans$coefsu, zcoefs = ans$zcoefs, pvals.beta = ans$pz, se.beta = ans$sez, gcv = ans$gcv, gcvu = ans$gcvu, xnms = xnms, znms = znms, zid = zid, vals = vals, zid1 = zid1, zid2 = zid2, ynm = ynm, decrs = dcs, tms = mt, is_param = is_param, is_fac = is_fac)
  rslt$call <- cl
  class(rslt) <- "wps"
  return (rslt) 
}

#####
wps.fit = function(x1t, x2t, y, zmat = NULL, w = NULL, pen = 0, pnt = FALSE, cpar = 1.5, decrs = c(FALSE, FALSE), delta = NULL, kts = NULL) {
	n = length(y)
	if (length(x1t) != n) {
		stop ("Error: length of x1 must equal length of y")
	}
	if (length(x2t) != n) {
		stop ("Error: length of x2 must equal length of y")
	}
# new
	k1 = kts[[1]]
	m1 = length(k1)
	k2 = kts[[2]]
	m2 = length(k2) 
	if (decrs[1]) {
		x1 = -x1t
		#if (!is.null(k1)) { 
			k1 = -k1[m1:1] 
		#}
	} else {x1 = x1t}
	if (decrs[2]) {
		x2 = -x2t
		#if (!is.null(k2)) {
			k2 = -k2[m2:1]
		#}
	} else {x2 = x2t}
	one = 1:n*0 + 1
	if (is.null(zmat)) {
		zmat = matrix(one, ncol = 1)
	} else {
		if (dim(zmat)[1] != n) {
			stop ("Error: # rows of zmat must equal length of y")
		}
		zproj = zmat %*% solve(crossprod(zmat), t(zmat))
		onep = one - zproj %*% one
# if the one vector is not in the space of zmat, then include the one vector in zmat 
		if (sum(onep^2) > 1e-12) {
			zmat = cbind(one, zmat)
		}
	}
	p = dim(zmat)[2]
	dimnames(zmat)[[2]] = NULL
# we already get warp.delta; the 1st col of xmat is one
	xmat0 = delta
#print (all(xmat0[, 1] == 1))
	if (dim(zmat)[2] > 1) {
		xmat = cbind(zmat, xmat0[, 2:(m1 + m2 - 1 + (m1 - 1) * (m2 - 1))])
	} else {xmat = xmat0}
# constraint matrix
	amat = matrix(0, nrow = 2*m1*m2 - m1 - m2, ncol = m1*m2)
	amat[1, 2] = 1
	amat[2, m1 + 1] = 1
	irow = 2
	for (i2 in 2:m2) {
		irow = irow+1
		amat[irow, 2] = 1; amat[irow, m1 + m2 + i2 - 2] = 1
	}
	for (i1 in 2:m1) {
		irow = irow + 1
		amat[irow, m1 + 1] = 1; amat[irow, m1 + m2 + (i1 - 2) * (m2 - 1)] = 1
	}
	for (i1 in 2:(m1 - 1)) {
		irow = irow + 1
		amat[irow, i1] = -1; amat[irow, i1 + 1] = 1
	}
	for (i2 in 2:(m2 - 1)) {
		irow = irow + 1
		amat[irow, m1 + i2 - 1] = -1
		amat[irow, m1 + i2] = 1
	}
	for (i1 in 2:(m1 - 1)) {
		for (i2 in 2:m2) {
			irow = irow + 1
			amat[irow, i1] = -1
			amat[irow, i1 + 1] = 1
			amat[irow, m1 + m2 - 2 + (i1 - 2) * (m2 - 1) + i2] = -1
			amat[irow, m1 + m2 - 2 + (i1 - 1) * (m2 - 1) + i2] = 1
		}
	}
	for (i2 in 2:(m2 - 1)) {
		for (i1 in 2:m1) {
			irow = irow + 1
			amat[irow, m1 + i2 - 1] = -1
			amat[irow, m1 + i2] = 1
			amat[irow, m1 + m2 - 2 + (i1 - 2) * (m2 - 1) + i2] = -1
			amat[irow, m1 + m2 - 1 + (i1 - 2) * (m2 - 1) + i2] = 1
		}
	}
	amat0 = amat
	if (dim(zmat)[2] > 1) {
		amatz =  matrix(0, nrow = 2 * m1 * m2 - m1 - m2, ncol = dim(zmat)[2])
		amat = cbind(amatz, amat[, 2:(m1 * m2)])
	}
# penalty matrix
	kn1 = k1; kn2 = k2
	dmat = matrix(0, nrow = 2 * (m1 * m2 - m1 - m2), ncol = m1 * m2)
# by row
	irow = 1
	dmat[irow, 2] = -2; dmat[irow, 3] = 1
	for (i in 2:(m1 - 2)) {
		irow = irow + 1
		dmat[irow, i] = 1; dmat[irow, i + 1] = -2; dmat[irow, i + 2] = 1
	}
	for (ik2 in 2:m2) {
		irow = irow + 1
		dmat[irow, 2] = -2; dmat[irow, 3] = 1
		dmat[irow, m1 + m2 - 1 + ik2 - 1] = -2
		dmat[irow, m1 + 2 * m2 + ik2 - 3] = 1
		for (ik1 in 2:(m1 - 2)) {
			irow = irow + 1
			dmat[irow, ik1] = 1; dmat[irow, ik1 + 1] = -2; dmat[irow, ik1 + 2] = 1
			dmat[irow, m1 + m2 - 1 + (m2 - 1) * (ik1 - 2) + ik2 - 1] = 1
			dmat[irow, m1 + m2 - 1 + (m2 - 1) * (ik1 - 1) + ik2 - 1] = -2
			dmat[irow, m1 + m2 - 1 + (m2 - 1) * ik1 + ik2 - 1] = 1
		}
	}
# by col
	irow = irow + 1
	dmat[irow, m1 + 1] = -2; dmat[irow, m1 + 2] = 1
	for (i in 2:(m2 - 2)) {
		irow = irow + 1
		dmat[irow, m1 + i - 1] = 1; dmat[irow, m1 + i] = -2; dmat[irow, m1 + i + 1] = 1
	}
	for (ik1 in 2:m1) {
		irow = irow + 1
		dmat[irow, m1 + 1] = -2; dmat[irow, m1 + 2] = 1
		dmat[irow, m1 + m2 - 1 + (m2 - 1) * (ik1 - 2) + 1] = -2
		dmat[irow, m1 + m2 - 1 + (m2 - 1) * (ik1 - 2) + 2] = 1
		for (ik2 in 2:(m2 - 2)) {
			irow = irow + 1
			dmat[irow, m1 + ik2 - 1] = 1; dmat[irow, m1 + ik2] = -2; dmat[irow, m1 + ik2 + 1] = 1
			dmat[irow, m1 + m2 - 1 + (m2 - 1) * (ik1 - 2) + ik2 - 1] = 1
			dmat[irow, m1 + m2 - 1 + (m2 - 1) * (ik1 - 2) + ik2] = -2
			dmat[irow, m1 + m2 - 1 + (m2 - 1) * (ik1 - 2) + ik2 + 1] = 1
		}
	}
	dmat0 = dmat
	if (dim(zmat)[2] > 1) {
		dmatz = matrix(0, nrow = dim(dmat)[1], ncol = dim(zmat)[2])
		dmat = cbind(dmatz, dmat[, 2:(m1 * m2)])
	}
	k1 = kn1; k2 = kn2
#}
# weight
	if (is.null(w)) {
		xw0 = delta
		yw = y
		xw = xmat
		zw = zmat
	} else {
		if (min(w) > 1e-8) {
			yw = y * sqrt(w)
			xw = xmat
			zw = zmat
			xw0 = xmat0
			for (i in 1:n) {
				xw[i, ] = xmat[i, ] * sqrt(w[i])
				xw0[i, ] = xmat0[i, ] * sqrt(w[i])
			}
			for (i in 1:n) {
				zw[i, ] = zmat[i, ] * sqrt(w[i])
			}
		} else {
			xw0 = xmat0; xw = xmat; yw = y; zw = zmat
		}
	}
# transform to cone projection
	ps = max(1e-6, pen)
	qmat = t(xw) %*% xw + ps * t(dmat) %*% dmat	
	qmat0 = t(xw0) %*% xw0 + ps * t(dmat0) %*% dmat0
	umat = chol(qmat)
	#iumat = diag(ncol(umat))
        #uinv = backsolve(umat, iumat)
	#umat0 = chol(qmat0)
	uinv = solve(umat)
	#uinv0 = solve(umat0)
	atil = amat %*% uinv
	#atil0 = amat0 %*% uinv0
	cvec = t(uinv) %*% t(xw) %*% yw
	ans = coneA(cvec, atil)
	phihat = ans$thetahat
	ahat = uinv %*% phihat
	#cvec = crossprod(xw, yw)
	#cvec = t(xw) %*% yw	
#bv = rep(0, nrow(amat))
	#ans = qprog(qmat, cvec, amat, 1:nrow(amat)*0, msg = FALSE)
	#ahat = ans$thetahat
	muhat = xw %*% ahat	
# get trace of "proj" matrix
	sm = 1e-10
	nz = 1:dim(amat)[1] < sm	
	nz[amat %*% ahat > sm] = TRUE
	if (sum(nz) < dim(amat0)[1]) {
		gamj = -t(amat0[!nz, ])
		if (length(gamj) == dim(xw0)[2]) {
			gamj = matrix(gamj, ncol = 1)
		}
		aq = qr(gamj)
		if (aq$rank >= 1) {
			if (aq$rank == 1) {
				gamj = matrix(qr.Q(aq)[, 1:aq$rank, drop = FALSE], ncol = 1)
			} else {
				gamj = qr.Q(aq)[, 1:aq$rank, drop = FALSE]
			}	
			#pa = gamj %*% solve(crossprod(gamj), t(gamj))		
			pa = gamj %*% solve(t(gamj) %*% gamj) %*% t(gamj)	
			dj1 = dim(gamj)[1]; dj2 = dim(gamj)[2]
			imat = matrix(0, nrow = dj1, ncol = dj1) 
			for (i in 1:dj1) {imat[i, i] = 1}
			uvecs = matrix(rnorm(dj1 * (dj1 - dj2)), nrow = dj1)
			wmatt = (imat - pa) %*% uvecs
			aqt = qr(wmatt); wmat = qr.Q(aqt)
			b0mat = xw0 %*% wmat
		}
	} else {
		b0mat = xw0
		#inmat = matrix(0, nrow = dim(xmat0)[2], ncol = dim(xmat0)[2])
		#for (i in 1:dim(xmat0)[2]) {inmat[i, i] = 1}
		inmat = diag(dim(xmat0)[2])
		wmat = inmat		
	}
	p0mat = xmat0 %*% wmat %*% solve(t(wmat) %*% qmat0 %*% wmat) %*% t(wmat) %*% t(xmat0)
	#p0mat = xmat0 %*% wmat %*% solve(crossprod(wmat, qmat0 %*% wmat), t(wmat) %*% t(xmat0))	
	edf = sum(diag(p0mat)) + p - 1
	inmat = matrix(0, nrow = n, ncol = n)
	diag(inmat) = 1	
	sse1 = sum((yw - xw %*% ahat)^2)
#new:
zcoefs = ahat[1:p]
sse0 = sum((yw - zw %*% zcoefs)^2)
#new use edf instead if (n - cpar * edf) < 0
if ((n - cpar * edf) <= 0) {
	sig2hat = sse1 / edf
} else {
	sig2hat = sse1 / (n - cpar * edf)
}
	if (p >= 1) {
		pm = one %*% solve(crossprod(one)) %*% t(one)
		covmat = solve(t(zmat) %*% (inmat - p0mat + pm) %*% zmat)	
		sez = sqrt(diag(covmat) * sig2hat)	
		tz = zcoefs / sez
#new use edf instead if (n - cpar * edf) < 0
		if ((n - cpar * edf) <= 0) {
			pz = 2 * (1 - pt(abs(tz), edf))
			warning ('Effective degrees of freedom is close to the number of observations! Inference about parametric covariates is not reliable!')
#print ('Check pz!')
		} else {
			pz = 2 * (1 - pt(abs(tz), n - cpar * edf))
		}
	}
# get unconstrained penalized estimator
	#prmatu = xw %*% solve(qmat, t(xw))
	prmatu = xw %*% solve(qmat) %*% t(xw)
	muhatu = prmatu %*% yw
	sseu = sum((yw - muhatu)^2)
	edfu = sum(diag(prmatu))
	ans = new.env()
	ans$k1 = k1
	ans$k2 = k2
	ans$muhat = muhat
	#ans$muplot = mupl
	ans$muhatu = muhatu
	#ans$muplotu = muplu
	ans$gcv = sse1 / (1 - edf / n)^2
	ans$gcvu = sseu / (1 - edfu/n)^2
	#ans$ssr = sse1
	ans$sse1 = sse1
	ans$sse0 = sse0
	ans$edf = edf
	#ans$nz = amat %*% ahat
	ans$edfu = edfu
	ans$coefs = ahat
	ans$coefsu = solve(qmat) %*% t(xw) %*% yw
#include coef for one vector
	#ans$zcoefs = ahat[1:p]
	ans$zcoefs = zcoefs
	ans$zmat = zmat 	
	ans$sighat = sig2hat
	ans$delta = xmat
	ans$pen = pen
	#ans$gcvus = gcvus
	#ans$lambdas_pen = lambdas_pen
	#if (p >= 2) {
	if (p >= 1) {
		ans$sez = sez
		ans$pz = pz
		ans$covmat = covmat
		ans$tz = tz
	} else {ans$sez = NULL; ans$pz = NULL; ans$covmat = NULL; ans$tz = NULL}
#print (sig2hat)
	return (ans)
}

####################################################################
#four monotonicity functions for warped-plane fit
####################################################################
ii <- function(x1, x2, numknots = c(0, 0), knots = list(k1 = 0, k2 = 0), space = c("E", "E"))
{
    cl <- match.call()
    pars1 <- match.call()[2]
    pars2 <- match.call()[3]
    xm <- cbind(x1, x2)
    attr(xm, "name") <- c(deparse(pars1$x1), deparse(pars2$x2))
    attr(xm, "shape") <- "wps_ii"
    attr(xm, "numknots") <- numknots
    attr(xm, "knots") <- knots
    attr(xm, "space") <- space
    attr(xm, "decreasing") <- c(FALSE, FALSE)
    return (xm)
}


di <- function(x1, x2, numknots = c(0, 0), knots = list(k1 = 0, k2 = 0), space = c("E", "E"))
{
    cl <- match.call()
    pars1 <- match.call()[2]
    pars2 <- match.call()[3]
    xm <- cbind(x1, x2)
    attr(xm, "name") <- c(deparse(pars1$x1), deparse(pars2$x2))
    attr(xm, "shape") <- "wps_di"
    attr(xm, "numknots") <- numknots
    attr(xm, "knots") <- knots
    attr(xm, "space") <- space
    attr(xm, "decreasing") <- c(TRUE, FALSE)
    return (xm)
}


dd <- function(x1, x2, numknots = c(0, 0), knots = list(k1 = 0, k2 = 0), space = c("E", "E"))
{
    cl <- match.call()
    pars1 <- match.call()[2]
    pars2 <- match.call()[3]
    xm <- cbind(x1, x2)
    attr(xm, "name") <- c(deparse(pars1$x1), deparse(pars2$x2))
    attr(xm, "shape") <- "wps_dd"
    attr(xm, "numknots") <- numknots
    attr(xm, "knots") <- knots
    attr(xm, "space") <- space
    attr(xm, "decreasing") <- c(TRUE, TRUE)
    return (xm)
}


###############################################################
#makedelta_wps: make delta to a pair of warped-plane variables#
###############################################################
makedelta_wps <- function(x1t, x2t, m1_0 = 0, m2_0 = 0, k1 = 0, k2 = 0, space = c("E", "E"), decreasing = c(FALSE, FALSE))
{
# x1 and x2 no need to sort
# if decreasing  (all calculations done for doubly-increasing case)
	n = length(x1t)
	if (decreasing[1]) {
#print (k1 != 0)
		x1 = -x1t
		#if (!is.null(k1)) {
		if (k1 != 0) { 
			m1 = length(k1); k1 = -k1[m1:1] 
		}
	} else {x1 = x1t}
	if (decreasing[2]) {
		x2 = -x2t
		#if (!is.null(k2)) {
		if (k2 != 0) {
			m2 = length(k2); k2 = -k2[m2:1]
		}
	} else {x2 = x2t}
# determine whether to use default knots or user-defined knots
# check if k1 and k2 > 1 elems
	make1 = FALSE   
	#if (is.null(k1)) {
	#if (k1 == 0) {
	if (all(k1 == 0)) {
		make1 = TRUE
	} else if (min(k1) > min(x1) | max(k1) < max(x1)) {
		warning ('Predictor should be within the range of the knots! User-defined knots not used!')
		make1 = TRUE
# new:	at least 4 kts?
	}# else if (length(k1) < 4) {
	#	make1 = TRUE
	#} else if (m1_0 >= 4) {
	#	make1 = TRUE
	#}
	make2 = FALSE
	#if (is.null(k2)) {
	#if (k2 == 0) {
	if (all(k2 == 0)) {
		make2 = TRUE
	} else if (min(k2) > min(x2) | max(k2) < max(x2)) {
		warning ('Predictor should be within the range of the knots! User-defined knots not used!')
		make2 = TRUE
	}# else if (length(k2) < 4) {
	#	make2 = TRUE
	#} else if (m2_0 >= 4) {
	#	make2 = TRUE
	#}
# add quantile part in it
	if (make1) {
# new:
#print (m1_0 == 12L)? !is.integer(m1_0) not work
		if (m1_0  < 4 | round(m1_0, 0) != m1_0) {
		#if (m1_0  < 4) {
			warning ('At least four knots should be used! Number of knots is re-defined!')
			m1 = 2 * round(n^(1/6)) + 4 
		} else {m1 = m1_0}
		if (space[1] == "Q") {
			k1 = quantile(unique(x1), probs = seq(0, 1, length = m1), names = FALSE)
		}
		if (space[1] == "E") {
			k1 = 0:(m1 - 1) / (m1 - 1) * (max(x1) - min(x1)) + min(x1)
		} 
		#k1 = 0:(m1 - 1) / (m1 - 1) * (max(x1) - min(x1)) + min(x1)
	} else { m1 = length(k1) }
	if (make2) {
#new:
		if (m2_0 < 4 | round(m2_0, 0) != m2_0) {
			warning ('At least four knots should be used! Number of knots is re-defined!')
			m2 = 2 * round(n^(1/6)) + 4 
		} else {m2 = m2_0}
		if (space[2] == "Q") {
			k2 = quantile(unique(x2), probs = seq(0, 1, length = m2), names = FALSE)
		}
		if (space[2] == "E") {
			k2 = 0:(m2 - 1) / (m2 - 1) * (max(x2) - min(x2)) + min(x2)
		} 
		#k2 = 0:(m2 - 1) / (m2 - 1) * (max(x2) - min(x2)) + min(x2)	
	} else { m2 = length(k2) }
## check to see if empty knot intervals
	rm1 = k1[m1]; rm2 = k2[m2]
	keep = 1:m1 > 0
	for (i in 2:m1) {
		if (sum(x1 >= k1[i-1] & x1 < k1[i]) == 0) {
			keep[i] = FALSE
		}
	}
	k1 = k1[keep]; m1 = length(k1); k1[m1] = rm1
	keep = 1:m2 > 0
	for (i in 2:m2) {
		if (sum(x2 >= k2[i-1] & x2 < k2[i]) == 0) {
			keep[i] = FALSE
		}
	}
	k2 = k2[keep]; m2 = length(k2); k2[m2] = rm2
# make the basis functions
	b1 = matrix(0, nrow = n, ncol = m1)
	b2 = matrix(0, nrow = n, ncol = m2)
	for (i in 2:(m1 - 1)) {
		i1 = x1 >= k1[i-1] & x1 <= k1[i]
		b1[i1, i] = (x1[i1] - k1[i-1]) / (k1[i] - k1[i-1])
		i2 = x1 > k1[i] & x1 <= k1[i+1]
		b1[i2, i] = (k1[i+1] - x1[i2]) / (k1[i+1] - k1[i])
	}
	i1 = x1 >= k1[1] & x1 <= k1[2]
	b1[i1, 1] = (k1[2] - x1[i1]) / (k1[2] - k1[1])
	i2 = x1 > k1[m1-1] & x1 <= k1[m1]
	b1[i2, m1] = (x1[i2] - k1[m1-1]) / (k1[m1] - k1[m1-1])

	for (i in 2:(m2 - 1)) {
		i1 = x2 >= k2[i-1] & x2 <= k2[i]
		b2[i1, i] = (x2[i1] - k2[i-1]) / (k2[i] - k2[i-1])
		i2 = x2 > k2[i] & x2 <= k2[i+1]
		b2[i2, i] = (k2[i+1] - x2[i2]) / (k2[i+1] - k2[i])
	}
	i1 = x2 >= k2[1] & x2 <= k2[2]
	b2[i1, 1] = (k2[2] - x2[i1]) / (k2[2] - k2[1])
	i2 = x2 > k2[m2-1] & x2 <= k2[m2]
	b2[i2, m2] = (x2[i2] - k2[m2-1]) / (k2[m2] - k2[m2-1])
## design matrix
	xmat0 = matrix(nrow = n, ncol = m1 + m2 - 1 + (m1 - 1) * (m2 - 1))
	xmat0[ ,1] = 1:n*0 + 1
	xmat0[ ,2:m1] = b1[ ,2:m1]
	xmat0[ ,(m1 + 1):(m1 + m2 - 1)] = b2[ ,2:m2]
	for (i in 1:(m1 - 1)) {
		xmat0[ ,(m1 + m2 + (i - 1) * (m2 - 1)):(m1 + m2 - 1 + i * (m2 - 1))] = b1[ ,i + 1] * b2[ ,2:m2]
	}
	#if (dim(zmat)[2] > 1) {
	#	xmat = cbind(zmat, xmat0[ ,2:(m1 + m2 - 1 + (m1 - 1) * (m2 - 1))])
	#} else {xmat = xmat0}
	#bmat = xmat0[ ,2:(m1 + m2 - 1 + (m1 - 1) * (m2 - 1))]
# ignore the zmat for now
	xmat = xmat0
# columns of delta are edges, different from other make_delta in cgam	
	delta = xmat
	ans = new.env()
	ans$delta = delta
	ans$k1 = k1
	ans$k2 = k2
	return (ans)	
	#attr(delta, "shape") = "warp"
	#delta
}













