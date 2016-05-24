#' Functions copied and modified from survival package
#' @description The distcomp package makes use of code from the survival package,
#' with the permission of the original authors. This includes R code as well as C code.
#' That is, the underlying Cox model code is derived from that in the R survival package.
#' The original copyrights are retained for these files and the notices preserved.
#' However, these are for internal use and future implementations may change how we use them.
#' In order to avoid confusion and any name collision, the names of these functions have
#' been modified to include a prefix "dc".
#'
#' @import survival
#' @importFrom stats terms
#' @importFrom stats model.extract
#' @importFrom stats model.matrix
#' @importFrom stats model.offset
#' @importFrom stats model.weights
#' @importFrom stats .getXlevels
#' @useDynLib distcomp, .registration = TRUE
#' @keywords internal
#tt <- function(x) x
dccoxph <- function(formula, data, weights, subset, na.action,
                    init, control, ties= c("efron", "breslow", "exact"),
                    singular.ok =TRUE, robust=FALSE,
                    model=FALSE, x=FALSE, y=TRUE,  tt, method=ties, ...) {

    ties <- match.arg(ties)
    Call <- match.call()

                                        # create a call to model.frame() that contains the formula (required)
                                        #  and any other of the relevant optional arguments
                                        # then evaluate it in the proper frame
    indx <- match(c("formula", "data", "weights", "subset", "na.action"),
                  names(Call), nomatch=0)
    if (indx[1] ==0) stop("A formula argument is required")
    temp <- Call[c(1,indx)]  # only keep the arguments we wanted
    temp[[1]] <- as.name('model.frame')  # change the function called

    special <- c("strata", "cluster", "tt")
    temp$formula <- if(missing(data)) terms(formula, special)
    else              terms(formula, special, data=data)
                                        # Make "tt" visible for coxph formulas, without making it visible elsewhere
    if (!is.null(attr(temp$formula, "specials")$tt)) {
        coxenv <- new.env(parent= environment(formula))
        assign("tt", function(x) x, env=coxenv)
        environment(temp$formula) <- coxenv
    }

    mf <- eval(temp, parent.frame())
    if (nrow(mf) ==0) stop("No (non-missing) observations")
    Terms <- terms(mf)


    ## We want to pass any ... args to coxph.control, but not pass things
    ##  like "dats=mydata" where someone just made a typo.  The use of ...
    ##  is simply to allow things like "eps=1e6" with easier typing
    extraArgs <- list(...)
    if (length(extraArgs)) {
        controlargs <- names(formals(coxph.control)) #legal arg names
        indx <- pmatch(names(extraArgs), controlargs, nomatch=0L)
        if (any(indx==0L))
            stop(gettextf("Argument %s not matched", names(extraArgs)[indx==0L]),
                 domain = NA)
    }
    if (missing(control)) control <- coxph.control(...)

    ## Force iter.max = 0 to ensure only one iteration
    control$iter.max = 0L

    Y <- model.extract(mf, "response")
    if (!inherits(Y, "Surv")) stop("Response must be a survival object")
    type <- attr(Y, "type")
    if (type!='right' && type!='counting')
        stop(paste("Cox model doesn't support \"", type,
                   "\" survival data", sep=''))
    weights <- model.weights(mf)
    data.n <- nrow(Y)   #remember this before any time transforms

    cluster<- attr(Terms, "specials")$cluster
    if (length(cluster)) {
        robust <- TRUE  #flag to later compute a robust variance
        tempc <- untangle.specials(Terms, 'cluster', 1:10)
        ord <- attr(Terms, 'order')[tempc$terms]
        if (any(ord>1)) stop ("Cluster can not be used in an interaction")
        cluster <- strata(mf[,tempc$vars], shortlabel=TRUE)  #allow multiples
        dropterms <- tempc$terms  #we won't want this in the X matrix
        dropcon <- tempc$vars
                                        # Save away xlevels after removing cluster (we don't want to save upteen
                                        #  levels of that variable, which we will never need).
        xlevels <- .getXlevels(Terms[-tempc$terms], mf)
    } else {
        dropterms <- dropcons <- NULL
        if (missing(robust)) robust <- FALSE
        xlevels <- .getXlevels(Terms, mf)
    }

    strats <- attr(Terms, "specials")$strata
    if (length(strats)) {
        stemp <- untangle.specials(Terms, 'strata', 1)
        if (length(stemp$vars)==1) strata.keep <- mf[[stemp$vars]]
        else strata.keep <- strata(mf[,stemp$vars], shortlabel=TRUE)
        strats <- as.numeric(strata.keep)
    }

    timetrans <- attr(Terms, "specials")$tt
    if (missing(tt)) tt <- NULL
    if (length(timetrans)) {
        timetrans <- untangle.specials(Terms, 'tt')
        ntrans <- length(timetrans$terms)

        if (is.null(tt)) {
            tt <- function(x, time, riskset, weights){ #default to O'Brien's logit rank
                obrien <- function(x) {
                    r <- rank(x)
                    (r-.5)/(.5+length(r)-r)
                }
                unlist(tapply(x, riskset, obrien))
            }
        }
        if (is.function(tt)) tt <- list(tt)  #single function becomes a list

        if (is.list(tt)) {
            if (any(!sapply(tt, is.function)))
                stop("The tt argument must contain function or list of functions")
            if (length(tt) != ntrans) {
                if (length(tt) ==1) {
                    temp <- vector("list", ntrans)
                    for (i in 1:ntrans) temp[[i]] <- tt[[1]]
                    tt <- temp
                }
                else stop("Wrong length for tt argument")
            }
        } else stop("The tt argument must contain a function or list of functions")

        if (ncol(Y)==2) {
            if (length(strats)==0) {
                sorted <- order(-Y[,1], Y[,2])
                newstrat <- rep.int(0L, nrow(Y))
                newstrat[1] <- 1L
            } else {
                sorted <- order(strats, -Y[,1], Y[,2])
                                        #newstrat marks the first obs of each strata
                newstrat <-  as.integer(c(1, 1*(diff(strats[sorted])!=0)))
            }
            if (storage.mode(Y) != "double") storage.mode(Y) <- "double"
            counts <- .Call(Ccoxcount1, Y[sorted,],
                            as.integer(newstrat))
            tindex <- sorted[counts$index]
        }
        else {
            if (length(strats)==0) {
                sort.end  <- order(-Y[,2], Y[,3])
                sort.start<- order(-Y[,1])
                newstrat  <- c(1L, rep(0, nrow(Y) -1))
            } else {
                sort.end  <- order(strats, -Y[,2], Y[,3])
                sort.start<- order(strats, -Y[,1])
                newstrat  <- c(1L, as.integer(diff(strats[sort.end])!=0))
            }
            if (storage.mode(Y) != "double") storage.mode(Y) <- "double"
            counts <- .Call(Ccoxcount2, Y,
                            as.integer(sort.start -1L),
                            as.integer(sort.end -1L),
                            as.integer(newstrat))
            tindex <- counts$index
        }
        mf <- mf[tindex,]
        Y <- Surv(rep(counts$time, counts$nrisk), counts$status)
        type <- 'right'  # new Y is right censored, even if the old was (start, stop]
        strats <- rep(1:length(counts$nrisk), counts$nrisk)
        weights <- model.weights(mf)
        for (i in 1:ntrans)
            mf[[timetrans$var[i]]] <- (tt[[i]])(mf[[timetrans$var[i]]], Y[,1], strats,
                                                weights)
    }

    contrast.arg <- NULL  #due to shared code with model.matrix.coxph
    attr(Terms, "intercept") <- TRUE
    adrop <- 0  #levels of "assign" to be dropped; 0= intercept
    stemp <- untangle.specials(Terms, 'strata', 1)
    if (length(stemp$vars) > 0) {  #if there is a strata statement
        hasinteractions <- FALSE
        for (i in stemp$vars) {  #multiple strata terms are allowed
                                        # The factors att has one row for each variable in the frame, one
                                        #   col for each term in the model.  Pick rows for each strata
                                        #   var, and find if it participates in any interactions.
            if (any(attr(Terms, 'order')[attr(Terms, "factors")[i,] >0] >1))
                hasinteractions <- TRUE
        }
        if (!hasinteractions)
            dropterms <- c(dropterms, stemp$terms)
        else adrop <- c(0, match(stemp$var, colnames(attr(Terms, 'factors'))))
    }

    if (length(dropterms)) {
        temppred <- attr(terms, "predvars")
        Terms2 <- Terms[ -dropterms]
        if (!is.null(temppred)) {
                                        # subscripting a Terms object currently drops predvars, in error
            attr(Terms2, "predvars") <- temppred[-(1+dropterms)] # "Call" object
        }
        X <- model.matrix(Terms2, mf, constrasts=contrast.arg)
                                        # we want to number the terms wrt the original model matrix
                                        # Do not forget the intercept, which will be a zero
        renumber <- match(colnames(attr(Terms2, "factors")),
                          colnames(attr(Terms,  "factors")))
        attr(X, "assign") <- c(0, renumber)[1+attr(X, "assign")]
    }
    else X <- model.matrix(Terms, mf, contrasts=contrast.arg)

                                        # drop the intercept after the fact, and also drop strata if necessary
    Xatt <- attributes(X)
    xdrop <- Xatt$assign %in% adrop  #columns to drop (always the intercept)
    X <- X[, !xdrop, drop=FALSE]
    attr(X, "assign") <- Xatt$assign[!xdrop]
                                        #if (any(adrop>0)) attr(X, "contrasts") <- Xatt$contrasts[-adrop]
                                        #else attr(X, "contrasts") <- Xatt$contrasts
    attr(X, "contrasts") <- Xatt$contrasts
    offset <- model.offset(mf)
    if (is.null(offset) | all(offset==0)) offset <- rep(0., nrow(mf))

    assign <- attrassign(X, Terms)
    contr.save <- attr(X, "contrasts")
    if (missing(init)) init <- NULL
    pterms <- sapply(mf, inherits, 'coxph.penalty')
    if (any(pterms)) {
        pattr <- lapply(mf[pterms], attributes)
        pname <- names(pterms)[pterms]
                                        #
                                        # Check the order of any penalty terms
        ord <- attr(Terms, "order")[match(pname, attr(Terms, 'term.labels'))]
        if (any(ord>1)) stop ('Penalty terms cannot be in an interaction')
        pcols <- assign[match(pname, names(assign))]
        list(p = ncol(X),
             fitter = function(beta) {
                 fit <- do.call(getFromNamespace("coxpenal.fit", "survival"),
                                list(X, Y, strats, offset, init=beta,
                                     control,
                                     weights=weights, method=method,
                                     row.names(mf), pcols, pattr, assign))
                 fit
             })
    } else {
        if( method=="breslow" || method =="efron") {
            if (type== 'right')  fitter <- get("dccoxph.fit")
            else                 fitter <- get("agreg.fit")
        } else if (method=='exact') {
            if (type== "right")  fitter <- get("coxexact.fit")
            else  fitter <- get("agexact.fit")
        }  else stop(paste ("Unknown method", method))

        list(p = ncol(X),
             fitter = function(beta) {
                 fit <- fitter(X, Y, strats, offset, init=beta, control, weights=weights,
                               method=method, row.names(mf))
                 fit
             })
    }
}

#' @rdname dccoxph
dccoxph.fit <- function(x, y, strata, offset, init, control,
			weights, method, rownames)
    {
        n <-  nrow(y)
        if (is.matrix(x)) nvar <- ncol(x)
        else {
            if (length(x)==0) nvar <-0
            else nvar <-1
	}
        time <- y[,1]
        status <- y[,2]

                                        # Sort the data (or rather, get a list of sorted indices)
        if (length(strata)==0) {
            sorted <- order(time)
            strata <- NULL
            newstrat <- as.integer(rep(0,n))
	}
        else {
            sorted <- order(strata, time)
            strata <- strata[sorted]
            newstrat <- as.integer(c(1*(diff(as.numeric(strata))!=0), 1))
	}
        if (missing(offset) || is.null(offset)) offset <- rep(0,n)
        if (missing(weights)|| is.null(weights))weights<- rep(1,n)
        else {
            if (any(weights<=0)) stop("Invalid weights, must be >0")
            weights <- weights[sorted]
	}
        stime <- as.double(time[sorted])
        sstat <- as.integer(status[sorted])

        if (nvar==0) {
                                        # A special case: Null model.
                                        #  (This is why I need the rownames arg- can't use x' names)
                                        # Set things up for 0 iterations on a dummy variable
            x <- as.matrix(rep(1.0, n))
            nullmodel <- TRUE
            nvar <- 1
            init <- 0
            maxiter <- 0
	}
        else {
            nullmodel <- FALSE
            maxiter <- control$iter.max
            if (!missing(init) && length(init)>0) {
                if (length(init) != nvar) stop("Wrong length for inital values")
	    }
            else init <- rep(0,nvar)
	}

        storage.mode(weights) <- storage.mode(init) <- "double"
        coxfit <- .Call(Ccoxfit6,
                        as.integer(maxiter),
                        stime,
                        sstat,
                        x[sorted,],
                        as.double(offset[sorted]),
                        weights,
                        newstrat,
                        as.integer(method=="efron"),
                        as.double(control$eps),
                        as.double(control$toler.chol),
                        as.vector(init),
                        as.integer(1))  # internally rescale

        if (nullmodel) {
            score <- exp(offset[sorted])
            coxres <- .C(Ccoxmart, as.integer(n),
                         as.integer(method=='efron'),
                         stime,
                         sstat,
                         newstrat,
                         as.double(score),
                         as.double(weights),
                         resid=double(n))

            resid <- double(n)
            resid[sorted] <- coxres$resid
            names(resid) <- rownames

            list( loglik = coxfit$loglik[1],
                 linear.predictors = offset,
                 residuals = resid,
                 method= c('coxph.null', 'coxph') )
        }
        else {
            var <- matrix(coxfit$imat,nvar,nvar)
            coef <- coxfit$coef
            if (coxfit$flag < nvar) which.sing <- diag(var)==0
            else which.sing <- rep(FALSE,nvar)

            infs <- abs(coxfit$u %*% var)
            if (maxiter >1) {
                if (coxfit$flag == 1000)
                    warning("Ran out of iterations and did not converge")
                else {
                    infs <- ((infs > control$eps) &
                                 infs > control$toler.inf*abs(coef))
                    if (any(infs))
                        warning(paste("Loglik converged before variable ",
                                      paste((1:nvar)[infs],collapse=","),
                                      "; beta may be infinite. "))
		}
	    }

            names(coef) <- dimnames(x)[[2]]
            lp <- c(x %*% coef) + offset - sum(coef*coxfit$means)
            score <- exp(lp[sorted])
            coxres <- .C(Ccoxmart, as.integer(n),
                         as.integer(method=='efron'),
                         stime,
                         sstat,
                         newstrat,
                         as.double(score),
                         as.double(weights),
                         resid=double(n))
            resid <- double(n)
            resid[sorted] <- coxres$resid
            names(resid) <- rownames
            if (maxiter > 0) coef[which.sing] <- NA  #leave it be if iter=0 is set

            concordance <- survConcordance.fit(Surv(stime, sstat), lp[sorted],
                                               strata, weights)
            list(coefficients  = coef,
                 var    = var,
                 loglik = coxfit$loglik,
                 score  = coxfit$sctest,
                 iter   = coxfit$iter,
                 linear.predictors = as.vector(lp),
                 residuals = resid,
                 means = coxfit$means,
                 imatCopy = coxfit$imatCopy,
                 concordance=concordance,
                 gradient = coxfit$u,
                 method='coxph')
	}
    }
