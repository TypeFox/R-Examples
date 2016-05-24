#ALG: itree modifies this substantially from the rpart functions,
# which is called 'rpart'.

# SCCS  @(#)rpart.s	1.35 07/05/01
#
#  The recursive partitioning function, for S
#
# ALG 4/4/2012: updated to accept penalty and
# impscale as arguments.  impscale is whether we
# scale the improvement by root or parent.
# Make the default settings correspond to
# original rpart.

itree <- function(formula, data, weights, subset,
		   na.action=na.itree, method, penalty="none",
		   model=FALSE, x=FALSE, y=TRUE,
		   parms, control, cost, ...)
{
    call <- match.call()

    #check if we've been passed a dataframe in the model argument
    if (is.data.frame(model)) {
		m <- model
		model <- FALSE
	}
    else {
    	#no? create the model ourselves, set everything to
    	#null and fill in later. Added penalty & impscale 4/4/2012
		m <- match.call(expand.dots=FALSE)
		m$model <- m$method <- m$control<- NULL
		m$x <- m$y <- m$parms <- m$... <- NULL
		m$penalty <- NULL
		m$cost <- NULL
		m$na.action <- na.action
		m[[1L]] <- as.name("model.frame")
		m <- eval(m, parent.frame())
	}

    Terms <- attr(m, "terms")

    #check if the model has been passed interactions.
	#If yes, stop and complain.
    if(any(attr(Terms, "order") > 1L)){
		stop("Trees cannot handle interaction terms")
	}


    Y <- model.extract(m, "response")
    wt <- model.extract(m, "weights")

    #If no weights have been passed, set them all to 1.0
    if(length(wt)==0L){
    	wt <- rep(1.0, nrow(m))
    }

    offset <- attr(Terms, "offset")

    #design matrix is X
    X <- itree.matrix(m)
    nobs <- nrow(X)
    nvar <- ncol(X)

	# Figure out the method = min(SS),purity, gini, etc.
    if (missing(method)) {
	if (is.factor(Y) || is.character(Y))      method <- 'class'
        else if (inherits(Y, "Surv"))   method <- 'exp'
	else if (is.matrix(Y)) method<- 'poisson'
	else                   method<- 'anova'
	}

    if (is.list(method)) {
	# User written split methods
		mlist <- method
		method <- 'user'

		if (missing(parms)) init <- mlist$init(Y, offset, wt=wt)
		else                init <- mlist$init(Y, offset, parms, wt)

		method.int <- 4L     #the 4th in func_table_objective.h

	        ## assign this to avoid garbage collection
	        keep <- itreecallback(mlist, nobs, init)
    }
    # pre-defined split
    else {

    	#ALG 9/19/2012:
    	# for method 'purity' or 'extremes' we need to figure out if its class/quantitative
    	#so we know which purity or extremes functions to call.
		if(method %in% c("purity","extremes")){
			#check that it's not survival data.
			if (inherits(Y, "Surv")) {stop("This method isn't defined for Survival data.")}

			# 'y' categorical --> reset method
			if (is.factor(Y) || is.character(Y)){
				#method is class_purity or class_extremes
				method <- paste("class",method,sep="_")
			}
			else{
				#regression_purity or regression_extremes
				method <- paste("regression",method,sep="_")
			}
		}
		method.int <- pmatch(method, c("anova", "poisson", "class", "exp",
							"regression_extremes","regression_purity",
							"class_extremes","class_purity"))

		if (is.na(method.int)) stop("Invalid method")
		method <- c("anova", "poisson", "class", "exp",
					"regression_extremes","regression_purity","class_extremes","class_purity")[method.int]

		#ALG 6/19/2013: temporarilty turn of exp and poisson in this version
		# since they're not checked yet.
		if(method %in% c("exp","poisson")){
			stop("That method is not currently available in itree.")
		}
		if (method.int==4L) method.int <- 2L  #exp goes to poisson.

		if (missing(parms))
		  init <- (get(paste("itree", method, sep='.')))(Y, offset, ,wt)
		else
		  init <- (get(paste("itree", method, sep='.')))(Y, offset, parms, wt)

	        ns <- asNamespace("itree")
	        if(!is.null(init$print)) environment(init$print) <- ns
	        if(!is.null(init$summary)) environment(init$summary) <- ns
	        if(!is.null(init$text)) environment(init$text) <- ns
    }

	# ALG 4/3/2012: FIGURE OUT THE PENALTY TYPE
	# default is none.
	penalty.int <- pmatch(penalty, c("none","newvar","ema"))
	if (is.na(penalty.int)) stop("Invalid penalty")
	penalty <- c("none","newvar","ema")[penalty.int]

    Y <- init$y

    xlevels <- attr(X, "column.levels")
    cats <- rep(0,ncol(X))
    if(!is.null(xlevels)) {
	cats[match(names(xlevels), dimnames(X)[[2]])] <-
		  unlist(lapply(xlevels, length))
	}

    # We want to pass any ... args to rpart.control, but not pass things
    #  like "dats=mydata" where someone just made a typo.  The use of ...
    #  is just to allow things like "cp=.05" with easier typing
    extraArgs <- list(...)
    indx <- 0
    if (length(extraArgs)) {
		controlargs <- names(formals(itree.control))  #legal arg names
		indx <- match(names(extraArgs), controlargs, nomatch=0)

		if (any(indx==0))
            stop("Argument ", names(extraArgs)[indx==0], "not matched")
	}

    controls <- itree.control(...)
    if (!missing(control)) {
    #not missing... check that method and cp don't collide
    	controls[names(control)] <- control
    }

	#ALG 4/28/2012: choose default cp if none passed,
	#otherwise make sure cp & method combination is ok.
	cp_passed <- (sum(indx==12)>0)+(!is.null(controls$cp))
	if(!cp_passed){
		controls$cp <- .01
		if(method.int %in% c(5,6,7,8)){
			controls$cp <- 0
		}
	}
	else{
		#check that what was passed doesn't conflict with the method.
		if( (method.int %in% c(5,6,7,8)) && controls$cp!=0)
			stop("This method is only defined for cp=0.")
	}

    xval <- controls$xval
    if (is.null(xval) || (length(xval)==1L && xval==0) || method=='user') {
		xgroups <-0
		xval <- 0
	}
    else if (length(xval)==1L) {
	# make random groups
        xgroups <- sample(rep(1:xval, length=nobs), nobs, replace=FALSE)
	}
    else if (length(xval) == nobs) {
	xgroups <- xval
	xval <- length(unique(xgroups))
	}
    else {
	# Check to see if observations were removed due to missing
	if (!is.null(attr(m, 'na.action'))) {
	    # if na.rpart was used, then na.action will be a vector
	    temp <- as.integer(attr(m, 'na.action'))
	    xval <- xval[-temp]
	    if (length(xval) == nobs) {
		xgroups <- xval
		xval <- length(unique(xgroups))
		}
	    else stop("Wrong length for xval")
	    }
	else stop("Wrong length for xval")
	}

    #
    # Incorporate costs
    #
    # ALG 3/3/2012: this is where we enforce or throw an error if there
    # are costs passed in addition to asking for interpretable trees
    if (missing(cost)) cost <- rep(1.0, nvar)
    else {
		if (length(cost) != nvar){
			stop("Cost vector is the wrong length")
		}
		if (any(cost <=0)){
			 stop("Cost vector must be positive")
		}

		#ALG 7/25/2012. no cost for methods 5 & 6, 7, 8
		#no cost if penalty
		if( method.int %in% c(5,6,7,8)){
			stop("This method is not defined for variable costs")
		}
		if( penalty.int != 1){
			stop("Interpretability penalties are not defined for variable costs")
		}
	}

	# ALG 4/19/2012 if impscale is not specified, choose default
	# depending on the penalty & method
	# default is in a matrix with method along row and penalty along column
	# 1 <==> scale by parent, 2 <==> scale by root.
	default.impscale <- matrix(1,nrow=8,ncol=3)  	# Everything is scale by parent for now.
	if(controls$impscale==3){
		controls$impscale <- default.impscale[method.int, penalty.int]
	}

    #
    # Have s_to_rp consider ordered categories as continuous
    #  A right-hand side variable that is a matrix forms a special case
    # for the code.
    #
    # alg 2/11/2012: interp params go in 'controls'
    tfun <- function(x) {
	if (is.matrix(x)) rep(is.ordered(x), ncol(x))
	else is.ordered(x)
	}
    isord <- unlist(lapply(m[attr(Terms, 'term.labels')], tfun))
    rpfit <- .C(C_s_to_rp,
		    n = as.integer(nobs),
		    nvarx = as.integer(nvar),
		    ncat = as.integer(cats* !isord),
		    method= as.integer(method.int),
		    penalty=as.integer(penalty.int),
		    as.double(unlist(controls)),
		    parms = as.double(unlist(init$parms)),
		    as.integer(xval),
		    as.integer(xgroups),
		    as.double(t(init$y)),
		    as.double(X),
		    as.integer(!is.finite(X)), # R lets Infs through
		    error = character(1),
		    wt = as.double(wt),
		    as.integer(init$numy),
		    as.double(cost),
		    NAOK=TRUE)
    if (rpfit$n == -1)  stop(rpfit$error)

    # rpfit$newX[1:n] contains the final sorted order of the observations
    nodes <- rpfit$n          # total number of nodes
    nsplit<- rpfit$nvarx      # total number of splits, primary and surrogate
    numcp <- rpfit$method     # number of lines in cp table
    ncat  <- rpfit$ncat[1]    #total number of categorical splits
    numresp<- init$numresp    # length of the response vector

    if (nsplit == 0) xval <- 0
    cpcol <- if (xval>0 && nsplit>0) 5L else 3L
    if (ncat==0) catmat <- 0
    else         catmat <- matrix(integer(1), ncat, max(cats))

    rp    <- .C(C_s_to_rp2,
		       as.integer(nobs),
		       as.integer(nsplit),
		       as.integer(nodes),
		       as.integer(ncat),
		       as.integer(cats *!isord),
		       as.integer(max(cats)),
		       as.integer(xval),
		       which = integer(nobs),
		       cptable = matrix(double(numcp*cpcol), nrow=cpcol),
		       dsplit =  matrix(double(1),  nsplit,3),
		       isplit =  matrix(integer(1), nsplit,3),
		       csplit =  catmat,
		       dnode  =  matrix(double(1),  nodes, 3+numresp),
		       inode  =  matrix(integer(1), nodes, 6))
    tname <- c("<leaf>", dimnames(X)[[2]])
    #print(rp$dnode)

    if (cpcol==3) temp <- c("CP", "nsplit", "rel error")
    else          temp <- c("CP", "nsplit", "rel error", "xerror", "xstd")
    dimnames(rp$cptable) <- list(temp, 1L:numcp)

    # R change for empty-vector calculations.
    dn1 <- if(nsplit == 0L) character(0L) else tname[rp$isplit[,1L]+1L]
    splits<- matrix(c(rp$isplit[,2L:3L], rp$dsplit), ncol=5L,
                    dimnames = list(dn1,
                    c("count", "ncat", "improve", "index", "adj")))
    index <- rp$inode[,2]  #points to the first split for each node

    # Now, make ordered categories look like categories again (a printout
    #  choice)
    nadd <- sum(isord[rp$isplit[,1L]])
    if (nadd >0) {
	newc <- matrix(integer(1), nadd, max(cats))
	cvar <- rp$isplit[,1L]
	indx <- isord[cvar]		     # vector of T/F
	cdir <- splits[indx,2L]               # which direction splits went
	ccut <- floor(splits[indx,4L])        # cut point
	splits[indx,2L] <- cats[cvar[indx]]   #Now, # of categories instead
	splits[indx,4L] <- ncat + 1L:nadd      # rows to contain the splits

	# Next 4 lines can be done without a loop, but become indecipherable
	for (i in 1L:nadd) {
	    newc[i, 1L:(cats[(cvar[indx])[i]])] <- -1*as.integer(cdir[i])
	    newc[i, 1L:ccut[i]] <- as.integer(cdir[i])
	    }
	if (ncat==0) catmat <- newc
	else         catmat <- rbind(rp$csplit, newc)
	ncat <- ncat + nadd
	}
    else catmat <- rp$csplit

    if (nsplit==0) {  #tree with no splits
	frame <- data.frame(row.names=1,
			    var=  "<leaf>",
			    n =   rp$inode[,5L],
			    wt=   rp$dnode[,3L],
			    dev=  rp$dnode[,1L],
			    yval= rp$dnode[,4L],
			    complexity=rp$dnode[,2L],
			    ncompete  = pmax(0L, rp$inode[,3L] - 1L),
			    nsurrogate=rp$inode[,4L])
	}
    else {
	temp <- ifelse(index==0, 1, index)
	svar <- ifelse(index==0, 0, rp$isplit[temp,1L]) #var number
	frame <- data.frame(row.names=rp$inode[,1L],
			    var=  factor(svar, 0:ncol(X), tname),
			    n =   rp$inode[,5L],
			    wt=   rp$dnode[,3L],
			    dev=  rp$dnode[,1L],
			    yval= rp$dnode[,4L],
			    complexity=rp$dnode[,2L],
			    ncompete  = pmax(0L, rp$inode[,3L] - 1L),
			    nsurrogate=rp$inode[,4L])
	}
	#ALG: 9/24/2012: adjusted this to include all classification-based methods.
    if (method.int %in% c(3L,7L,8L)) {
        numclass <- init$numresp -1L
        # Create the class probability vector from the class counts, and
        #   add it to the results
        # The "pmax" one line down is for the case of a factor y which has
        #   no one at all in one of its classes.  Both the prior and the
        #   count will be zero, which led to a 0/0.
        temp <- rp$dnode[,-(1L:4L), drop = FALSE] %*% diag(init$parms$prior*
                                           sum(init$counts)/pmax(1,init$counts))
        yprob <- temp /rowSums(temp)   #necessary with altered priors
        yval2 <- matrix(rp$dnode[, -(1L:3L)], ncol=numclass+1)
		yval2 <- as.data.frame(cbind(yval2[,-1],yprob)) #first col is just a repeat of yval.
		colnames(yval2) <- c(paste("wt.class",1:numclass,sep=""),paste("wt.frac.class",1:numclass,sep=""))
        nodewt <- frame$wt/frame$wt[1]
		#frame$yval2 <- cbind(yval2, yprob)
		#frame$nodewt <- nodewt
		frame <- cbind(frame,yval2,nodewt)
	}
    else if (init$numresp >1L) frame$yval2 <- rp$dnode[,-(1L:3L), drop = FALSE]

    if (is.null(init$summary))
	    stop("Initialization routine is missing the summary function")
    if (is.null(init$print))
	    functions <- list(summary=init$summary)
    else    functions <- list(summary=init$summary, print=init$print)
    if (!is.null(init$text)) functions <- c(functions, list(text=init$text))
    if (method=='user')	functions <- c(functions, mlist)

    where <- rp$which
    names(where) <- row.names(m)

    if (nsplit ==0L) {  # no 'splits' component
	ans <- list(frame = frame,
		    where = where,
		    call=call, terms=Terms,
		    cptable =  t(rp$cptable),
		    method = method,
		    penalty = penalty,
		    parms  = init$parms,
		    control= controls,
		    functions= functions)
	}
    else {
	ans <- list(frame = frame,
		    where = where,
		    call=call, terms=Terms,
		    cptable =  t(rp$cptable),
		    splits = splits,
		    method = method,
		    penalty = penalty,
		    parms  = init$parms,
		    control= controls,
		    functions= functions)

		# remove cp table if it's a one-sided method.
		if(method.int %in% c(5,6,7,8)){
			ans$cptable <- NULL
		}
	}
    if (ncat>0) ans$csplit <- catmat + 2L
    if (model) {
	ans$model <- m
	if (missing(y)) y <- FALSE
	}
    if (y) ans$y <- Y
    if (x) {
	ans$x <- X
	ans$wt<- wt
	}
    ans$ordered <- isord
    if(!is.null(attr(m, "na.action")))
	    ans$na.action <- attr(m, "na.action")
    if (!is.null(xlevels)) attr(ans, 'xlevels') <- xlevels

    if (method.int %in% c(3L,7L,8L)){
    	attr(ans, "ylevels") <- init$ylevels
    }

# # if (length(xgroups)) ans$xgroups <- xgroups
    class(ans) <- "itree"
    ans
 }
