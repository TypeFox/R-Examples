#### $Id: qbsks.R,v 1.22 2009/02/24 13:47:41 maechler Exp $

qbsks2 <-
    function(x,y,w,pw, knots,nknots, degree, Tlambda,
             constraint, ptConstr, maxiter, trace,
             nrq,nl1, neqc, tau, select.lambda,
             ks, do.select, knots.add, repeat.delete.add, ic, print.mesg,
             give.pseudo.x = TRUE,
             rq.tol = 1e-8, tol.kn = 1e-6, tol.0res = 1e-6, print.warn, nk.start)
{
    ##=########################################################################
    ##
    ## Compute B-spline coefficients for quantile B-spline with stepwise knots
    ## selection or with fixed knots (REGRESSION SPLINE), using
    ## Koenker and Ng (2005)  `Inequality Constrained Quantile Regression,"
    ## Sankhya, The Indian Journal of Statistics, 67, 418-440.
    ##
    ##=########################################################################

    smll.log <- 50*floor(log(.Machine$double.xmin)/50) # heavily rounding down
    finite.log <- function(x) {
        r <- log(x)
        if(is.na(r) || r > -Inf) r else smll.log # = -750 for IEEE arithmetic
    }
    n <- nrq
    ## xo <- x[order(x)]  <-- TODO ?
    logn <- log(n)
    f.IC <- switch(ic,
		   "AIC" = 2,
		   "BIC" =, "SIC" = logn,
		   stop("in qbsks2(): invalid 'ic' = ", ic, call. = FALSE))
    stopifnot(nk.start >= 2)
    constraint.orig <- constraint

    n.gr.sm <- with(ptConstr, n.greater + n.smaller)

    if(do.select) { ##  perform first step : knots selection
	if(length(Tlambda) > 1)
	    stop("you cannot select knots for more than one lambda simultaneously")
        if(print.mesg) cat("qbsks2():\n Performing general knot selection ...\n")#4
        Tic <- Tifl <- double(nknots-1)
        for(i in (nk.start-1):(nknots-1)) {
            Tknots <- knots[seq(1,nknots, len = i+1)]
            n.Tknts <- length(Tknots)
            if(n.Tknts == 2 && degree == 1 && constraint %in% c("convex", "concave"))
                ## guard against trying convex fit when only 2 knots are used
                constraint <- "none"
            dim.o <- getdim2(degree, n.Tknts, constraint)
            ks <- dim.o$ks
            Tnvar <- dim.o$nvar
            niqc <- dim.o$n.iqc + n.gr.sm
	    rqss <- drqssbc2(x,y,w, pw, Tknots, degree, Tlambda,
			     constraint, ptConstr, maxiter, trace=trace-1,
			     nrq,nl1,neqc,niqc, Tnvar,
			     tau=tau, select.lambda=select.lambda,
                             give.pseudo.x = give.pseudo.x,
			     rq.tol = rq.tol, tol.0res = tol.0res,
                             print.warn=print.warn)
            constraint <- constraint.orig
            Tic[i] <- finite.log(rqss$fidel) -logn + (i-1+ks)*f.IC / n
            Tifl[i] <- rqss$ifl
        }

        ## "FIXME": allow  "+ 1 S.E." selection instead of ``simple arg_min''
        Tic.min <- min(Tic)
        Tifl.final <- Tifl[Tic == Tic.min]
        nknots.min <- min((1:(nknots-1))[Tic == Tic.min])
        if(nknots.min == 1 && Tifl.final != 1) {
            ##
            ## when the chosen nknots=2, guard against anomaly of ifl=5 when
            ## constraint=='periodic', or ifl=2 when the chosen model is infeasible.
            ##
            Tic.min <- min(Tic[2:length(Tic)])
            Tifl.final <- Tifl[Tic == Tic.min]
            nknots.min <- min((1:(nknots-1))[Tic == Tic.min])
        }
        if(Tifl.final %in% c(2,3,4))
            return(list(ifl = Tifl.final))

	warnUP <- function(nknots, ic)
	    ## warn(5,nknots,ic)
	    cat("\n WARNING! Since the number of ",nknots," knots selected by ",
		ic," reached the\n",
		"  upper bound during general knot selection, you might want to rerun\n",
		"  cobs with a larger number of knots. \n")

        up.bound.reached <- (nknots.min + 1 == nknots)
	if(up.bound.reached && print.warn)
	    warnUP(nknots, ic)

        knots <- knots[seq(1,nknots, len = nknots.min+1)]
        names(knots) <- NULL
        ##
        ## perform knots deletion
        ##
        repeat.a.d <- TRUE
        Tic.global.min <- Tic.min
        while(repeat.a.d) {
            delete <- TRUE
            if(print.mesg) cat("\n Deleting unnecessary knots ...\n") # 5
            while(delete && nknots.min > 1) {
                n.Tknts <- length(knots)
                Tic1 <- rep.int(0, n.Tknts-2)
                n.Tknts.1 <- n.Tknts - 1
                if(n.Tknts.1 == 2 && degree == 1 && constraint %in% c("convex", "concave"))
                    ## guard against convex fit when only 2 knots are used
                    constraint <- "none"
                dim.o <- getdim2(degree,n.Tknts.1, constraint)
                ks <- dim.o$ks
                Tnvar <- dim.o$nvar
                niqc <- dim.o$n.iqc + n.gr.sm
                for(i in 2:(n.Tknts-1)) {
                    Tknots <- knots[-i]
		    rqss <- drqssbc2(x, y, w, pw, Tknots, degree, Tlambda, constraint,
				     ptConstr, maxiter, trace=trace-1,
				     nrq, nl1, neqc, niqc, Tnvar,
				     tau=tau, select.lambda=select.lambda,
				     give.pseudo.x = give.pseudo.x,
				     rq.tol = rq.tol, tol.0res = tol.0res,
				     print.warn = print.warn)
                    constraint <- constraint.orig
                    Tic1[i-1] <- finite.log(rqss$fidel)-logn+(n.Tknts.1-2+ks)*f.IC / n
                }
                Tic1.min <- min(Tic1)
                idx.del <- min((2:(n.Tknts-1))[Tic1 == Tic1.min])
                if((delete <- Tic1.min <= Tic.min)) {
                    Tic.min <- Tic1.min
                    if(print.mesg >= 3)
                        cat("\n A knot at ",signif(knots[idx.del]),
                            " is deleted.\n") # 6
                    knots <- knots[-idx.del]
                    nknots.min <- length(knots)-1
                }
            }## end while(delete...)
            if(print.mesg >= 2) cat("\n No more knot to be deleted.\n") # 7
            ##
            ## perform knots addition
            ##
            if(knots.add) {
                add <- TRUE
                n.Tknts <- length(knots)
                if(print.mesg) cat("\n Searching for missing knots ...\n") # 8
                while(add && n.Tknts < nknots) {
                    Tic2 <- double(n.Tknts-1)
                    knots.add.1 <- (knots[1:(n.Tknts-1)]+knots[2:n.Tknts])/2
                    n.Tknts.1 <- n.Tknts + 1
                    dim.o <- getdim2(degree,n.Tknts.1,constraint)
                    ks <- dim.o$ks
                    Tnvar <- dim.o$nvar
                    niqc <- dim.o$n.iqc + n.gr.sm
                    for(i in 1:(n.Tknts-1)) {
                        Tknots <- sort(c(knots,knots.add.1[i]))
                        if(length(unique(cut00(x, Tknots))) != n.Tknts)
                            Tic2[i] <- Tic.min+1
                        else {
			    rqss <-
				drqssbc2(x,y,w,pw,Tknots,degree, Tlambda,
					 constraint, ptConstr, maxiter, trace=trace-1,
					 nrq,nl1,neqc,niqc,Tnvar,
					 tau=tau, select.lambda=select.lambda,
					 give.pseudo.x = give.pseudo.x,
					 rq.tol = rq.tol, tol.0res = tol.0res,
					 print.warn = print.warn)

                            Tic2[i] <- finite.log(rqss$fidel) -logn +
                                (n.Tknts.1-2+ks)*f.IC / n
                        }
                    }
                    Tic2.min <- min(Tic2)
                    idx.add <- min((1:(n.Tknts-1))[Tic2 == Tic2.min])
                    if((add <- Tic2.min <= Tic.min)) {
                        Tic.min <- Tic2.min
                        knots <- sort(c(knots,knots.add.1[idx.add]))
                        if(print.mesg >= 2)
                            cat("\n A knot at ",signif(knots.add.1[idx.add]),
                                " is added.\n") # 9
                    }
                    n.Tknts <- length(knots)
                }## end while(add ..)
                if(print.mesg >= 2) cat("\n No more knot to be added.\n") # 10
                if(n.Tknts == nknots) up.bound.reached <- TRUE
            }## (knots.add)

            repeat.a.d <- repeat.delete.add && Tic.global.min > Tic.min
            if(repeat.a.d) {
                Tic.global.min <- Tic.min
            }
        }## end while(repeat.a.d)

	if(up.bound.reached && print.warn)
	    warnUP(nknots, ic)
        if(print.mesg >= 2) cat("\n Computing the final fit ...\n") # 11
    } ## end if(do.select)

    ##
    ## compute the B-spline coefficients for the full sample
    ##
    nknots <- length(knots)
    ## shift the first and last knot a tiny bit "outside":
    rk <- diff(range(knots))
    knots[1] <- knots[1] - tol.kn*rk
    knots[nknots] <- knots[nknots] + tol.kn*rk

    if(nknots == 2 && (constraint %in% c("convex", "concave")) &&
       degree == 1) { # guard against convex fit when only 2 knots are used
        if(print.warn)
	    cat(sprintf("WARNING: %s from \"%s\" to \"none\"\n",
			"only two knots; changing 'constraint'", constraint))
        constraint <- "none"
    }
    dim.o <- getdim2(degree, nknots, constraint)
    ks <- dim.o$ks
    Tnvar <- dim.o$nvar
    niqc <- dim.o$n.iqc + n.gr.sm
    rqss <- drqssbc2(x,y,w,pw, knots, degree,Tlambda,
		     constraint, ptConstr, maxiter, trace=trace-1,
		     nrq,nl1, neqc,niqc, Tnvar,
		     tau=tau, select.lambda=select.lambda,
		     give.pseudo.x = give.pseudo.x,
		     rq.tol = rq.tol, tol.0res = tol.0res,
		     print.warn = print.warn)
    ## constraint <- constraint.orig

    c(rqss[c("coef", "fidel", "ifl", "icyc", "pseudo.x")],
      list(k = nknots-2+ks, knots = knots, nknots = nknots,
	   nvar = Tnvar, lambda = Tlambda))

} ## end qbsks()
