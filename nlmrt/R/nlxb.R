nlxb <- function(formula, start, trace = FALSE, data = NULL, 
    lower = -Inf, upper = Inf, masked = NULL, control = list(), 
    ...) {
    # A simplified and hopefully robust alternative to finding
    # the nonlinear least squares minimizer that causes
    # 'formula' to give a minimal residual sum of squares.
    # 
    # nls.mn is particularly intended to allow for the
    # resolution of very ill-conditioned or else near
    # zero-residual problems for which the regular nls()
    # function is ill-suited. It may also be a useful
    # pre-processor for nls().
    # 
    # J C Nash 2012-3-4 nashjc _at_ uottawa.ca
    # 
    # formula looks like 'y~b1/(1+b2*exp(-b3*T))' start MUST be
    # a vector where all the elements are named: e.g.,
    # start=c(b1=200, b2=50, b3=0.3) trace -- TRUE for console
    # output data is a data frame containing data for variables
    # used in the formula that are NOT the parameters. This
    # data may also be defined in the parent frame i.e.,
    # 'global' to this function lower is a vector of lower
    # bounds upper is a vector of upper bounds masked is a
    # character vector of names of parameters that are fixed.
    # control is a list of control parameters. These are: ...
    # 
    # ... will need to contain data for other variables that
    # appear in the formula and are defined in a parent frame
    # (Not sure how needed??)
    # 
    # This variant uses a qr solution without forming the sum
    # of squares and cross products t(J)%*%J
    # 
    # Function to display SS and point
    showpoint <- function(SS, pnum) {
        pnames <- names(pnum)
        npar <- length(pnum)
        cat("lamda:", lamda, " SS=", SS, " at")
        for (i in 1:npar) {
            cat(" ", pnames[i], "=", pnum[i])
        }
        cat(" ", feval, "/", jeval)
        cat("\n")
    }
    ######### get data from data frame if exists
    ######### print(str(data))
    if (!is.null(data)) {
        for (dfn in names(data)) {
            cmd <- paste(dfn, "<-data$", dfn, "")
            eval(parse(text = cmd))
        }
    } else stop("'data' must be a list or an environment")
    # ensure params in vector
    pnames <- names(start)
    start <- as.numeric(start)
    names(start) <- pnames
    # bounds
    npar <- length(start)  # number of parameters
    if (length(lower) == 1) 
        lower <- rep(lower, npar)
    if (length(upper) == 1) 
        upper <- rep(upper, npar)
    # ?? more tests on bounds
    if (length(lower) != npar) 
        stop("Wrong length: lower")
    if (length(upper) != npar) 
        stop("Wrong length: upper")
    if (any(start < lower) || any(start > upper)) 
        stop("Infeasible start")
    if (trace) {
        cat("formula: ")
        print(formula)
        cat("lower:")
        print(lower)
        cat("upper:")
        print(upper)
    }
    # controls
    ctrl <- list(watch = FALSE, phi = 1, lamda = 1e-04, offset = 100, 
        laminc = 10, lamdec = 4, femax = 10000, jemax = 5000, rofftest = TRUE, 
        smallsstest = TRUE)
     ##   maxlamda <- 1e+60) ## dropped 130709
    epstol <- (.Machine$double.eps) * ctrl$offset
    ncontrol <- names(control)
    nctrl <- names(ctrl)
    for (onename in ncontrol) {
        if (!(onename %in% nctrl)) {
            if (trace) 
                cat("control ", onename, " is not in default set\n")
        }
        ctrl[onename] <- control[onename]
    }
    if (trace) 
        print(ctrl)
    phiroot <- sqrt(ctrl$phi)
    lamda <- ctrl$lamda
    offset <- ctrl$offset
    laminc <- ctrl$laminc
    lamdec <- ctrl$lamdec  # save typing
    watch <- ctrl$watch
    femax <- ctrl$femax
    jemax <- ctrl$jemax
    # First get all the variable names:
    vn <- all.vars(parse(text = formula))
    # Then see which ones are parameters (get their positions
    # in the set xx
    pnum <- start  # may simplify later??
    pnames <- names(pnum)
    bdmsk <- rep(1, npar)  # set all params free for now
    maskidx <- which(pnames %in% masked)  # NOTE: %in% not == or order gives trouble
    if (length(maskidx) > 0 && trace) {
        cat("The following parameters are masked:")
        print(pnames[maskidx])
    }
    bdmsk[maskidx] <- 0  # fixed parameters
    if (trace) {
        parpos <- match(pnames, vn)
        datvar <- vn[-parpos]  # NOT the parameters
        for (dvn in datvar) {
            cat("Data variable ", dvn, ":")
            print(eval(parse(text = dvn)))
        }
    }
    if (is.character(formula)) {
        es <- formula
    } else {
        tstr <- as.character(formula)  # note ordering of terms!
        es <- paste(tstr[[2]], "~", tstr[[3]], "")
    }
    # Now separate the sides
    parts <- strsplit(as.character(es), "~")[[1]]
    if (length(parts) != 2) 
        stop("Model expression is incorrect!")
    lhs <- parts[1]
    rhs <- parts[2]
    # And build the residual at the parameters
    resexp <- paste(rhs, "-", lhs, collapse = " ")
    for (i in 1:npar) {
        # put parameters in separate variables
        joe <- paste(pnames[[i]], "<-", pnum[[i]])
        eval(parse(text = joe))
    }
    gradexp <- deriv(parse(text = resexp), names(start))  # gradient expression
    resbest <- with(data, eval(parse(text = resexp)))
    ssbest <- crossprod(resbest)
    ss0 <- ssbest + 1.0 # a base value
    ssminval <- ssbest * epstol^4
    if (trace) cat("ssminval =",ssminval,"\n")
    feval <- 1
    pbest <- pnum
    feval <- 1  # number of function evaluations
    jeval <- 0  # number of Jacobian evaluations
    if (trace) {
        cat("Start:")
        showpoint(ssbest, pnum)
        if (watch) 
            tmp <- readline("Continue")
    }
    if (length(maskidx) == npar) 
        stop("All parameters are masked")  # Should we return?
    ssquares <- .Machine$double.xmax  # make it big
    newjac <- TRUE  # set control for computing Jacobian
    eqcount <- 0
    roffstop <- FALSE
    smallstop <- FALSE
    while ((! roffstop) && (eqcount < npar) 
             && (feval <= femax) && (jeval <= jemax)
             && (! smallstop) ) {
        if (newjac) 
            {
                bdmsk <- rep(1, npar)
                bdmsk[maskidx] <- 0
                bdmsk[which(pnum - lower < epstol * (abs(lower) + 
                  epstol))] <- -3
                bdmsk[which(upper - pnum < epstol * (abs(upper) + 
                  epstol))] <- -1
                if (trace && watch) {
                  cat("bdmsk:")
                  print(bdmsk)
                }
                J0 <- with(data, eval(gradexp))
                Jac <- attr(J0, "gradient")
                jeval <- jeval + 1  # count Jacobians
                if (any(is.na(Jac))) 
                  stop("NaN in Jacobian")
                JTJ <- crossprod(Jac)
                # tmp<-readline("cont.")
                gjty <- t(Jac) %*% resbest  # raw gradient
                for (i in 1:npar) {
                  bmi <- bdmsk[i]
                  if (bmi == 0) {
                    gjty[i] <- 0  # masked
                    Jac[, i] <- 0
                  }
                  if (bmi < 0) 
                    {
                      if ((2 + bmi) * gjty[i] > 0) {
                        # free parameter
                        bdmsk[i] <- 1
                        if (trace) 
                          cat("freeing parameter ", i, " now at ", 
                            pnum[i], "\n")
                      } else {
                        gjty[i] <- 0  # active bound
                        Jac[, i] <- 0
                        if (trace) 
                          cat("active bound ", i, " at ", pnum[i], 
                            "\n")
                      }
                    }  # bmi
                }  # end for loop
                JTJnew<-crossprod(Jac)
                if (npar == 1) dee <- diag(as.matrix(sqrt(diag(crossprod(Jac)))))
                else dee <- diag(sqrt(diag(crossprod(Jac))))  # to append to Jacobian
            }  # end newjac
        lamroot <- sqrt(lamda)
        JJ <- rbind(Jac, lamroot * dee, lamroot * phiroot * diag(npar))  # build the matrix
        JQR <- qr(JJ)  # ??try
        rplus <- c(resbest, rep(0, 2 * npar))
        roff <- max(abs(as.numeric(crossprod(qr.Q(JQR), rplus))))/ss0
        if (trace) cat("roff =", roff,"  converged = ",(roff <= sqrt(epstol)),"\n")
        if (ctrl$rofftest && (roff <= sqrt(epstol))) roffstop <- TRUE
#        tmp <- readline('cont')
        delta <- try(qr.coef(JQR, -rplus))  # Note appended rows of y)
        if (class(delta) == "try-error") {
            stop("qr.coef failed")
            if (lamda < 1000 * .Machine$double.eps) 
                lamda <- 1000 * .Machine$double.eps
            lamda <- laminc * lamda
            newjac <- FALSE  # increasing lamda -- don't re-evaluate
            if (trace) 
                cat(" Equation solve failure\n")
            feval <- feval + 1  # count as a function evaluation to force stop
        } else {
            # solution OK
            if (trace) {
              cat("delta:")
              print(delta)
              cat("gjty:")
              print(gjty)
            }
            gproj <- crossprod(delta, gjty)
            gangle <- gproj/sqrt(crossprod(gjty) * crossprod(delta))
            gangle <- 180 * acos(sign(gangle)*min(1, abs(gangle)))/pi
            if (trace) 
                cat("gradient projection = ", gproj, " g-delta-angle=", 
                  gangle, "\n")
            if (is.na(gproj)) stop("gproj is NA")
            if (is.na(gproj) || (gproj > 0)) {
                # uphill direction -- should NOT be possible
                if (lamda < 1000 * .Machine$double.eps) 
                  lamda <- 1000 * .Machine$double.eps
                lamda <- laminc * lamda
                newjac <- FALSE  # increasing lamda -- don't re-evaluate
                if (trace) 
                  cat(" Uphill search direction\n")
                ## if (lamda > maxlamda) 
                ##   eqcount <- npar  # force stop on big lamda
            } else {
                # downhill
                delta[maskidx] <- 0
                delta <- as.numeric(delta)
                if (trace && watch) {
                  cat("delta:")
                  print(delta)
                }
                step <- rep(1, npar)
                # for (i in 1:npar){ step[i]<-0 if (delta[i]>0)
                # step[i]<-(upper[i]-pbest[i])/delta[i] if (delta[i]<0)
                # step[i]<-(lower[i]-pbest[i])/delta[i] # positive }
                for (i in 1:npar) {
                  bd <- bdmsk[i]
                  da <- delta[i]
                  # if (trace) cat(i,' bdmsk=',bd,' delta=',da,'\n')
                  if (bd == 0 || ((bd == -3) && (da < 0)) || 
                    ((bd == -1) && (da > 0))) {
                    delta[i] <- 0
                  } else {
                    if (delta[i] > 0) 
                      step[i] <- (upper[i] - pbest[i])/delta[i]
                    if (delta[i] < 0) 
                      step[i] <- (lower[i] - pbest[i])/delta[i]  # positive
                  }
                }
                # if (trace && watch) { cat('step:') print(step) }
                stepsize <- min(1, step[which(delta != 0)])
                if (trace) 
                  cat("Stepsize=", stepsize, "\n")
                if (stepsize < .Machine$double.eps) {
                  if (lamda < 1000 * .Machine$double.eps) 
                    lamda <- 1000 * .Machine$double.eps
                  lamda <- laminc * lamda
                  newjac <- FALSE  # increasing lamda -- don't re-evaluate
                  if (trace) 
                    cat(" Stepsize too small\n")
                } else {
                  # continue
                  pnum <- pbest + stepsize * delta  # adjust (note POSITIVE here, but not in nlsmn0
                  names(pnum) <- pnames  # NOT inherited through %*% !!!
                  eqcount <- length(which((offset + pbest) == 
                    (offset + pnum)))
                  if (eqcount < npar) {
                    for (i in 1:npar) {
                      joe <- paste(pnames[[i]], "<-", pnum[[i]])
                      eval(parse(text = joe))
                    }
                    feval <- feval + 1  # count evaluations
                    resid <- with(data, eval(parse(text = resexp)))
                    ssquares <- as.numeric(crossprod(resid))
                    if (is.na(ssquares)) 
                      ssquares <- .Machine$double.xmax
                    if (ssquares >= ssbest) {
                      if (lamda < 1000 * .Machine$double.eps) 
                        lamda <- 1000 * .Machine$double.eps
                      lamda <- laminc * lamda
                      newjac <- FALSE  # increasing lamda -- don't re-evaluate
                      if (trace) 
                        showpoint(ssquares, pnum)
                    } else {
                      lamda <- lamdec * lamda/laminc
                      if (trace) {
                        cat("<<")
                        showpoint(ssquares, pnum)
                      }
                      ssbest <- ssquares
                      if (ctrl$smallsstest) { smallstop<-(ssbest <= ssminval) }
                      resbest <- resid
                      pbest <- pnum
                      newjac <- TRUE
                    }  # reduced sumsquares
                  } else {
                    # end if equcount
                    if (trace) 
                      cat("No parameter change\n")
                  }
                }  # end stepsize not too small
            }  # end downhill
        }  # solution OK
        if (watch) 
            tmp <- readline("Cycle")
    }  # end main while loop
    pnum <- as.vector(pnum)
    names(pnum) <- pnames
    result <- list(resid = resbest, jacobian = Jac, feval = feval, 
        jeval = jeval, coefficients = pnum, ssquares = ssbest, lower=lower, upper=upper, maskidx=maskidx)
    class(result) <- "nlmrt"
    result
}
