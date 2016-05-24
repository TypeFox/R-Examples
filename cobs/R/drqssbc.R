#### $Id: drqssbc.R,v 1.42 2011/04/25 15:48:34 maechler Exp $

drqssbc2 <-
    function(x,y, w = rep.int(1,n), pw, knots, degree, Tlambda, constraint,
             ptConstr, maxiter = 100, trace = 0,
             nrq = length(x), nl1, neqc, niqc, nvar, tau = 0.50,
             select.lambda = length(Tlambda) > 1, give.pseudo.x = FALSE,
             rq.tol = 1e-8 * sc.y, tol.0res = 1e-6, print.warn = TRUE,
             rq.print.warn = FALSE)
{
    ##=########################################################################
    ##
    ## Estimate the B-spline coefficients for quantile *smoothing* spline, using
    ##		Ng (1996)  `An Algorithm for Quantile Smoothing Splines',
    ##		Computational Statistics & Data Analysis, 22, 99-118.
    ##
    ##=########################################################################
    n <- length(x)
    stopifnot((nj0 <- length(Tlambda)) >= 1, is.numeric(Tlambda),
              is.numeric(x), is.numeric(y), is.numeric(w),
              length(y) == n, length(w) == n,
              is.numeric(tau), length(tau) == 1, 0 <= tau, tau <= 1,
              maxiter == round(maxiter), maxiter >= 1,
              is.logical(select.lambda), is.logical(print.warn),
              is.logical(rq.print.warn))

    if(nj0 > 1) {
	if(!select.lambda)
	    warning("no lambda selection but still length(Tlambda) > 1 -- not fully tested yet")
	## MM [FIXME]: for that case, we now return *all* lambda results;
	##	  but this (particularly "cobs"-methods) is completely untested

	Tlambda <- sort(Tlambda)# needed for selection (messages and plots)
    }

    ## Scale tolerance by the "scale" (y)
    sc.y <- mean(abs(y - mean(y)))
    if(sc.y < 1e-12 * (my <- median(abs(y)))) {
        if(print.warn)
            cat("Warning: 'y's are almost constant ==> sc.y := scale(y) := 1e-12 median|y|")
        sc.y <- 1e-12 * my
    }
    tol.0res.y <- tol.0res * sc.y
    ## lambda.cut <- 10/n

    ##* Note: if  nrq != length(x) , e.g., in case of sub.sampling+ fit full
    ##*  if(select.lambda < 0)
    ##*      n <- nrq

    ## Y values of pointwise constraints in ptConstr need
    ## to be sorted according to the X values (consistent with those (l1|loo).design2):
    Tequal   <- if(ptConstr$n.equal   > 0)  ptConstr$equal[order(ptConstr$equal [,2]),3] # else NULL
    Tsmaller <- if(ptConstr$n.smaller > 0) -ptConstr$smaller [order(ptConstr$smaller [,2]),3]
    Tgreater <- if(ptConstr$n.greater > 0)  ptConstr$greater [order(ptConstr$greater [,2]),3]
    Tgradien <- if(ptConstr$n.gradient> 0)  ptConstr$gradient[order(ptConstr$gradient[,2]),3]
    Y.ptConstr <- c(Tsmaller, Tgreater, Tequal, -1*Tequal, Tgradien, -1*Tgradien)

    ## double matrix -- to contain "sol"ution :
    solNms <- c("lambda", "icyc", "ifl", "fidel", "sum|res|_s", "k")
    sol1 <- matrix(0., length(solNms), nj0, dimnames= list(solNms, NULL))
    allCoef <- matrix(0., nvar, nj0)

    if(trace >= 2)
	summSparse <- function(X, namX = deparse(substitute(X))) {
	    force(namX)
	    stopifnot(length(d <- dim(X)) == 2)
	    k <- length(X@ra)
	    sprintf("%s %4d x %d (nz = %d =^= %5.2g%%)",
		    namX, d[1], d[2], k, k/prod(d))
	}

    rqCtrl <- list(maxiter=maxiter, warn.mesg = rq.print.warn, small= rq.tol)
    for(ilam in 1:nj0) { ## for each lambda in Tlambda[] : `` grid search ''
	XX <-
	    if(degree == 1)
		l1.design2(x, w, constraint, ptConstr, knots, pw, nrq = n,
			   nl1, neqc, niqc, nvar, Tlambda[ilam])
	    else
		loo.design2(x, w, constraint,ptConstr, knots, pw, nrq = n,
			    nl1, neqc, niqc, nvar, Tlambda[ilam])
	Xeq <- XX$Xeq
	fieq <- XX$fieq
	if(fieq)
	    Xieq <- XX$Xieq

        if(trace >= 2) {
            if(ilam == 1 || trace >= 3) {
                ## only for the first lambda: dim() is same!
                ch1 <- sprintf("%3s.design2(%s): -> ", c("l1","loo")[degree],
                               if(trace >= 3)
                               sprintf("lam= %10.4g",Tlambda[ilam]) else '')
                cat(ch1, summSparse(Xeq), "\n")
                if(fieq) {
                    cat(sprintf("%*s", nchar(ch1) - 1, " "),
                        summSparse(Xieq), "\n")
                }
            }
            else cat(".", if(ilam == nj0) "\n")
        }

	niqc1 <- {
	    if (degree == 1 || Tlambda[ilam] == 0)  0
	    else 2*(length(knots) - 1) }
	##U   if(any(ibig <- abs(Xeq$ra) > big)) { ## make sure that drqssbc won't overflow
	##U	Xeq$ra[ibig] <- Xeq$ra[ibig] / big^0.5
	##U	warning("re-scaling ", sum(ibig), "values of Lp-design `Xeq'")
	##U   }
	##U   if(any(ibig <- abs(Xieq$ra) > big)) { ## make sure that drqssbc won't overflow
	##U	Xieq$ra[ibig] <- Xieq$ra[ibig] / big^0.5
	##U	warning("re-scaling ", sum(ibig), "values of Lp-design `Xieq'")
	##U   }

        ## FIXME: (Tnobs, n0) do *not* depend on lambda
        ##         --> keep outside the for() loop!
	Tnobs <- nrow(Xeq) + (if(fieq) nrow(Xieq) else 0)
                                        # Tnobs = number of pseudo-observations
	n0 <- Tnobs -n - with(ptConstr,
			      2*(n.equal + n.gradient) + n.smaller + n.greater)
        ## There are 2 inequality constraints for each equality one ==> '2 * (...)'
	Y <- c(y*w, rep.int(0, n0), Y.ptConstr)
	##    Yeq <- Y[1:(nrq+nl1+neqc)]
	Yeq <- Y[1:(nrq+nl1)]
	##    if(fieq) Yieq <- Y[(nrq+nl1+neqc+1):Tnobs]
	if(fieq) {
	    Yieq <- Y[(nrq+nl1+1):Tnobs]
	    Yeq.. <- c(Yeq, Yieq)
	    Xeq.. <- rbind(Xeq, Xieq)
	}
	## niqc2 <- Tnobs - niqc1

	rhs <- (1-tau)*Xeq[1:nrq, ]
	rhs <- ## MM: would be nice to have colSums() for "matrix.csr"
	    colSums(as.matrix(if(nl1 > 0) rbind(rhs,
						.5*Xeq[(nrq+1):(nrq+nl1),])
			      else rhs))
	z0 <-
	    if(fieq)
		rq.fit.sfnc(Xeq,Yeq, Xieq,Yieq, tau=tau, rhs=rhs, control=rqCtrl)
	    else rq.fit.sfn(Xeq,Yeq,		tau=tau, rhs=rhs, control=rqCtrl)
        ## rq.fit.sfn[c] : both in ../../quantreg/R/sfn.R
        ##  these call .Fortran("srqfn[c]", ..) in
        ##    -> ../../quantreg/src/srqfn.c and ..../srqfnc.c
        ##  these both call chlfct() in ../../quantreg/src/chlfct.c   (was *.f)
        ##  is built on misc.        in ../../quantreg/src/cholesky.c (was *.f)
	if(any(is.na(z0$coef)))
	    stop("The combination of 'constraint' and 'pointwise' is causing problems\n",
		 "for the algorithm with the current knot selection,\n knots : ",
		 capture.output(str(knots, vec.len = 10)),"\n",
		 "Check the feasibility of your 'constraint' and 'pointwise'\n",
		 "or use a different knot selection.\n\n If you are using automatic ",
		 "knot selection, try increasing 'nk.start' from its default value of 2.")

	pseudo.resid <- (if(fieq) Yeq.. else Yeq) -
	    c(as.matrix((if(fieq) Xeq.. else Xeq) %*% as.matrix.csr(z0$coef)))
        if(any(is.na(pseudo.resid))) ## not really seen this case...
	    warning(sprintf("pseudo residuals [c(%s)] are NA; (n,Tnobs)= (%d,%d)",
			    paste(which(is.na(pseudo.resid)), collapse= ","), n,Tnobs))
	resid <- pseudo.resid[1:n]
	aRes <- abs(pseudo.resid)
	k.idx <- c(1:nrq, if((nn <- nrq+nl1+niqc1+1) <= Tnobs) nn:Tnobs)
###	k.idx <- c(1:nrq) #PN: this makes more sense than the previous def. of effective dimensionality
        ##was k.idx <- 1:Tnobs ; k.idx <- c(k.idx[k.idx <= nrq], k.idx[k.idx >  nrq+nl1+niqc1])
        ## or          (1:Tnobs <= nrq) | (1:Tnobs > nrq+nl1+niqc)
	sol1[,ilam] <- c(Tlambda[ilam],
			 z0$it,
			 z0$ierr + 1, # =: ifl, i.e., success <==> ifl == 1 <==> ierr == 0
			 sum((tau-(resid < 0))* resid), # =: fidel
			 sum(aRes[(nrq+1):(nrq+nl1)])/Tlambda[ilam],
			 ## Determination of the "effective dimensionality"  'k' :
			 ##PN: too stringent: sum(aRes[k.idx] < tol.0res * mean(aRes[k.idx]))
			 max(degree+1, sum(aRes[k.idx] < tol.0res.y))) # {a kludge ..}
        allCoef[,ilam] <- z0$coef

    }## end for(ilam ..)

    if(trace == 1) ## only output at all:
        cat(sprintf("fieq=%s -> Tnobs = %d, n0 = %d, |ptConstr| = %d\n",
                    format(fieq), Tnobs, n0, length(Y.ptConstr)))

    ## MM: return *sparse* pseudo.x
    ##     Further: this is X~ of the last, not the best lambda
    pseudo.x <- if(give.pseudo.x) (if(fieq) Xeq.. else Xeq) ## else NULL

    if(select.lambda) {
	##
	## search for optimal lambda
	##
	is.ifl.1 <- sol1["ifl",] == 1
        ##PN: Weed out lambdas that give rise to interpolating fit
	notInt <- sol1["fidel",] > tol.0res.y
	i.keep <- is.ifl.1 & notInt
	n.keep <- sum(i.keep)
        ## <NEW: PN 2006-08-24, 08-28>
	## k.diff.tol <- 1
	i.mask <- rep(FALSE,nj0)
	if(any(sol1["k",i.keep] > sqrt(n)) & min(sol1["k",i.keep]) <= sqrt(n)) {
            ##PN: Use k0 to guide the lower bound of the lambda grid; the lower bound
            ##	is set when k0 declines by more than k.diff.tol while lambda decreases among the grid of i.keep lambdas
            k.up <- min((1:n.keep)[sol1["k",i.keep] <= sqrt(n)])
            i.k.up <- (1:nj0)[i.keep][k.up]
            i.cut <- i.k.up - 1
            i.mask[1:i.cut] <- TRUE
            if(sol1["k",i.k.up] != sol1["k",nj0]) {
                ##PN: exclude situations when the horizontal fit through discrete y
                ## ends up with more zeros than two, hence, a non-monotonically
                ## decreasing "k" in the right tail
                i.keep <- i.keep & !i.mask
                n.keep <- sum(i.keep)
            }
	}
        ## </NEW>

	if(n.keep == 0)
	    stop("The problem is degenerate for the range of lambda specified.")

        if(any(!is.ifl.1)) { ## had some errors in rq.fit.sfn*()
	    sol.err <- t(sol1[, !is.ifl.1, drop = FALSE])
	    if(print.warn) { ## had some errors in rq.fit.sfn*()
		cat(gettextf("WARNING: Some lambdas had problems in %s():\n",
			     if(fieq) "rq.fit.sfnc" else "rq.fit.sfn"))
		print( sol.err )
	    }
        } else sol.err <- NULL

	solx <- sol1[, i.keep, drop = FALSE]
	sicA <- sol1["fidel",]/n * n^(sol1["k",] /(2*n))
        ## round, so the subsequent selection is "less wavery":
	sic. <- signif(sicA[i.keep], 6)
	min.s <- min(sic.)
	i.min <- (1:n.keep)[sic. == min.s]
	if(sic.[n.keep] <= min.s && min.s < sic.[1]) { # min. at rightmost, i.e. roughest
	    i.best <- max(i.min)
	    wrn.flag <- 1
	}
	else if(sic.[1] <= min.s && min.s < sic.[n.keep]) { # min. at leftmost, i.e.smoothest
	    i.best <- min(i.min)
	    wrn.flag <- 2
	}
	else {
	    i.best <- ceiling(median(i.min)) # take middle of all minimal values
	    if(min.s == sic.[1] && min.s == sic.[n.keep])
		wrn.flag <- 3
	    else  ## "good" case
		wrn.flag <- 4
	}
	##
        ## trap infeasible solution when lam < 0
        ##
        if(print.warn) {
	    wrn1 <- "\n WARNING!  Since the optimal lambda chosen by SIC "
	    wrnPlot <- paste("plot() the returned object",
			     "(which plots 'sic' against 'lambda')")
	    wrn3 <- "and possibly consider doing one of the following:"
	    wrn4 <- "(1) reduce 'lambda.lo', increase 'lambda.hi', increase 'lambda.length' or all of the above;"
	    ns <- "\n   "
	    switch(wrn.flag,
		   cat(wrn1, "reached the smoothest possible fit at `lambda.hi', you should",
		       wrnPlot, wrn3, wrn4,
		       "(2) decrease the number of knots.\n", sep = ns),
		   cat(wrn1, "corresponds to the roughest possible fit, you should",
		       wrnPlot, wrn3, wrn4,
		       ## typically "increase"; sometimes decreasing is appropriate:
		       "(2) modify the number of knots.\n", sep = ns),
		   cat(wrn1, "rests on a flat portion, you might", wrnPlot,
		       "to see if you want to reduce 'lambda.lo' and/or increase 'lambda.ho'\n",
		       sep = ns),
		   cat('', "The algorithm has converged.  You might", wrnPlot,
			"to see if you have found the global minimum of the information criterion",
			"so that you can determine if you need to adjust any or all of",
			"'lambda.lo', 'lambda.hi' and 'lambda.length' and refit the model.\n",
		       sep = ns))
            ## if(solx["lambda",,i.best] == (Tlambda[i.keep])[length(Tlambda[i.keep])]
            ## && solx["k",i.best] > max(2,niqc2))
        }
	list(coef = allCoef[, i.keep, drop = FALSE][, i.best],
	     fidel = solx["fidel",],
	     k = as.integer(sol1["k",notInt]), #PN was min(solx["k",i.best], length(knots)-2+degree+1),
	     kk= as.integer(solx["k",i.best]),
	     ifl = as.integer(sol1["ifl",notInt]),
	     icyc = as.integer(sol1["icyc",notInt]), nvar = nvar,
	     lambda = solx["lambda",i.best],
	     pp.lambda = sol1["lambda",notInt], sic = log(sicA[notInt]),
	     sol.err = sol.err, flag = wrn.flag,
	     pseudo.x = pseudo.x, i.mask = i.mask[notInt])
    }
    else { ## not 'select.lambda' ---
	l2 <- as.list(as.data.frame(t(sol1)))
	for(nn in c("k", "ifl", "icyc")) l2[[nn]] <- as.integer(l2[[nn]])
	c(list(coef = drop(allCoef), nvar = nvar, pseudo.x = pseudo.x),
	  l2)
    }
} ## drqssbc2
