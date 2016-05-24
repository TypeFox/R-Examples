##
## The main attributeble risk function
##
attribrisk <- function(formula, data, weights, subset, na.action,
                  varmethod = c('jackknife', 'bootstrap', 'none'),
                  conf= .95, baseline, k=20, control, 
                  model=FALSE, x=FALSE, y=FALSE, ...) {

    Call <- match.call()   ## save a copy of the call
    varmethod <- match.arg(varmethod)

    # create a call to model.frame() that contains the formula (required)
    #  and any other of the relevant optional arguments
    # then evaluate it in the proper frame
    indx <- match(c("formula", "data", "weights", "subset", "na.action"),
                  names(Call), nomatch=0) 
    if (indx[1] ==0) stop("A formula argument is required")
    temp <- Call[c(1,indx)]  # only keep the arguments we wanted
    temp[[1]] <- as.name('model.frame')  # change the function called
    
    special <- c("strata", "expos")
    temp$formula <- if(missing(data)) terms(formula, special)
                    else              terms(formula, special, data=data)

    # Make "expos" visible for attribRisk formulas without making it 
    #   visible globally
    riskenv <- new.env(parent= environment(formula))
    assign("expos", function(x) x, envir=riskenv)
    environment(temp$formula) <- riskenv
    
    mf <- eval(temp, parent.frame())  #evaluate the model.frame call
    Terms <- terms(mf)
    xlevels <- .getXlevels(Terms, mf)

    ## Extract relevant variables
    Y <- model.response(mf)
    if (is.null(Y)) stop("formula does not have a response variable")
    if (is.factor(Y)) {
        if (length(levels(Y)) > 2) 
            stop("response must have only 2 levels")
        Y <- as.numeric(Y) -1
    }
    if (!all(Y==0 | Y==1)) 
        stop ("Response must be a 0/1 variable or a class with two levels")
    n <- length(Y)
    if (n==0) stop("No observations in the data set")
        
    weights <- model.weights(mf)
    if (length(weights)==0) weights <- rep(1.0, n)
    offset <- model.offset(mf)
    if (length(offset)==0) offset <- rep(0.0, n)
    

    ## We want to pass any ... args to the control function, but not pass things
    ##  like "dats=mydata" where someone just made a typo.  The use of ...
    ##  is simply to allow things like "xval=10" with easier typing
    extraArgs <- list(...)
    if (length(extraArgs)) {
        controlargs <- names(formals(coxph.control)) #legal arg names
        indx <- pmatch(names(extraArgs), controlargs, nomatch=0L)
        if (any(indx==0L))
            stop(paste("Argument %s not matched", names(extraArgs)[indx==0L]))
    }
    if (missing(control)) control <- attribrisk.control(...)

        
    ## Find the exposure variables
    expos <- untangle.specials(Terms, 'expos', 1)
    if (length(expos$var) >0) {
        number.of.exposures <- length(expos$var)
        exposure.names <- expos$var
    }
    else stop("There must be at least one exposure variable")

    # Did they put exposure into an interaction?  Not allowed
    temp2 <- attr(Terms, 'factors')[temp$var,, drop=FALSE] 
    if (any(temp2==1 & rep(attr(Terms, "order"), each=nrow(temp2)) >1))
        stop("exposure variables cannot be in an interaction")
        
    match.id <- NULL     #matched analysis flag
    strats <- attr(Terms, "specials")$strata
    if (length(strats)) {
        temp <- untangle.specials(Terms, 'strata', 1)
        dropx <- temp$terms
        if (length(temp$vars)==1) strata.keep <- mf[[temp$vars]]
        else stop("Only 1 strata term allowed")
        match.id <- as.numeric(strata.keep)
        
        # subscripting the Terms object messes up predvars and assign attributes
        #  but we don't use them
        Terms2 <- Terms[-dropx]
        attr(Terms2, 'intercept') <- 1 # intercept required
    }
    else Terms2 <- Terms

    X <- model.matrix(Terms2, mf)
    
    #  The AR estimate is associated with some change: all the smokers into
    #    non-smokers, or all cholesterol levels to < 210, or ...
    #  If a baseline data set was provided, it gives new covariate values for
    #    each subject.  If not, exposure variables are converted to their
    #    baseline values: reference level for a factor, smallest value if they
    #    are a character, FALSE if logical, or 0 if numeric.
    #
    # Create data frame m2, which has all exposures replaced with reference
    #  values.  The user can supply a baseline data set, or we use the
    #  above numbers.  The key value for AR is
    #     exp((original X - modified X) %*% beta)
    #  see equation 20 of the technical report.
    if (missing(baseline)) { 
        m2 <- mf
        m2[expos$var] <- lapply(mf[expos$var], function(x) {
            if (is.character(x)) x <- factor(x)
            if (is.factor(x)) x <- factor(rep(levels(x)[1], n),
                                         levels=levels(x))
            else if (is.logical(x)) rep(FALSE, n)
            else rep(0, n)
        })
    }                  
    else {
        ## First, does the baseline data contain the exposure variable names?
        ##  If we check before trying to use a formula on it we can
        ##  create a nicer error message.
        ##
        # The baseline data only needs exposure variables
        tform <- formula(paste("~", paste(expos$vars, collapse='+')))
        environment(tform) <- riskenv  # make expos() visible to tform

        indx <- is.na(match(all.vars(tform), names(baseline)))
        if (any(indx)) {
            stop(paste("Variable(s)", all.vars(tform)[indx],
                   "not found in the baseline data set"))
        }

        if (missing(subset))
            mbase <- model.frame(tform, data=baseline, xlev=xlevels,
                                 na.action=function(x) x)
        else mbase <- model.frame(tform, data=baseline, xlev=xlevels,
                                  na.action=function(x) x, subset=subset)

        # Delete any rows from mbase that were deleted due to missing in
        #  the original data mf
        if (nrow(mbase) !=1 && !is.null(attr(mf, 'na.action'))) {
            droprows <- as.integer(attr(mf, 'na.action')) #rows to toss  
            if (any(droprows) > nrow(mbase)) 
                stop("data and baseline have different observation counts")
            mbase <- mbase[-droprows,]
        }
        if (nrow(mbase) != nrow(mf) && nrow(mbase) !=1)
            stop("data and baseline have different observation counts")
        if (any(is.na(mbase)))
            stop("baseline data set cannot contain missing values")
        
        m2 <- mf
        for (i in 1:number.of.exposures) {
            temp <- mf[[expos$vars[i]]]
            if (class(temp) == 'factor' || class(temp) == 'character') {
                indx <- match(mbase[[i]], temp)
                if (any(is.na(indx)))
                    stop("A reference value is not in the orginal data")
                else  mbase[[i]] <- as.factor(temp)[indx] #propogate the levels
            }
            else if (class(temp) != class(mbase[[i]]) &&
                     !(is.numeric(temp) && is.numeric(mbase[[i]]))) #int/double
              stop("a baseline variable's class does not match that in data")
          
            m2[[expos$vars[i]]] <- mbase[[i]]
        }
    }
 
    xbase <- X - model.matrix(Terms2, data=m2, 
                              contrasts.arg =attr(X, 'contrasts'))
 
    if (length(strats)) {
      X <- X[,-1, drop=F]  #remove the intercept colum, for conditional fits
      xbase <- xbase[,-1, drop=F]
      Y <- Surv(rep(1, nrow(X)), Y)  #make it a survival object
    }

    ##
    ## Do the fit
    ##
    fit <- attribrisk.fit(X, Y, weights, offset, match.id,  
                          xbase, fit=TRUE)
 
    # Set up to do the variance
    #
    # The "k" parameter can be a single number, in which case the data is
    #  divided into that many groups, or it can be a vector of length n
    #  containing the group numbers.   The latter for the rare person who
    #  wants to control which subjects in which group directly.
    # If there are strata (match) then the grouping needs to keep parts of
    #  a match together.  Cross validation and bootstrap pull sets of matches.
    # There is one other odd special case: this is one of the few routines
    #  where data with integer case weights is used, i.e., a weight of "8"
    #  means 8 individuals like this.  Detect this special case as all 
    #  weights >=1, all weights integer, and at least one weight >1.  In this
    #  case we do simple jackknife, ignoring "k".
    casewt <- (all(weights>=1 & weights==floor(weights)) && max(weights) >1)
    if (casewt) {
        if (!missing(k) && k<n)
            stop("caseweights + grouped cross-validation not allowed")
        if (length(strats)) xgrp <- match.id
            else xgrp <- 1:n  
    }
    else {
        if (varmethod=="bootstrap") k <- n  # effectively ignore k
        if (length(k) ==1) {
            if (k==0) k <- n  #convient shortcut
            if (k <2) 
                stop("jackknife error cannot be based on less than 2 groups")

            # If the data is stratified, don't assign subgroups that break
            #   a strata
            if (length(strats)) {
                nstrat <- max(match.id)
                if (k >= nstrat) xgrp <- match.id
                else {
                    temp <- sample(1:floor(k), nstrat, replace=TRUE)
                    xgrp <- temp[match.id]
                }
            }
            else {
                if (k >=n) xgrp <- 1:n
                else xgrp <- sample(1:floor(k), n, replace=TRUE)
            }
        }
        else if (length(k) >=n) {
            # missing values nuisance
            temp <- attr(mf, "na.action")
            if (!is.na(temp)) k <- k[-temp]  #remove missing
            if (length(k) != n) stop("wrong length for 'k'")
            xgrp <- as.numeric(as.factor(k)) #turn it into integer groups
            if (length(strats)) {
                # Ensure that the user didn't break any matched sets
                temp <- tapply(xgrp, match.id, length)
                if (any(temp>1)) 
                    stop("cross-validation sets are not allowed to split a match")
            }
        }
        else stop("wrong length for k")
    }

    # Compute a variance
    if (varmethod=="jackknife")
        fit$var <- attribrisk_jack(fit$attribrisk, match.id, X, Y, weights,
                                offset, xbase, xgroup=xgrp, casewt)
    else if (varmethod=="bootstrap") {
        temp <- attribrisk_boot(match.id, X, Y, weights, offset, xbase,
                                casewt, xgroup=xgrp, 
                                conf, nboot=control$nboot, bci=control$bootci)
        fit$var <- temp$var
        fit$boot <- temp$boot
        fit$boot.ci  <- temp$boot.ci
    }
    fit$conf <- conf
    fit$terms <- Terms
    if(model)
      fit$model <- mf
    if(x)
      fit$x <- X
    if(y)
      fit$y <- Y

    if(!is.null(attr(mf, "na.action")))
      fit$na.action <- attr(mf, "na.action")
    fit$call <- Call

    class(fit) <- c("attribrisk", if (length(strats)) "coxph" else "glm")
    return(fit)
}
