mincrit <-function (obs, param, sumstats, obspar = NULL, abcmethod = abc, 
    crit = nn.ent, sumsubs = 1:ncol(sumstats), limit = length(sumsubs), 
    do.only = NULL, verbose = TRUE, do.crit = TRUE, do.err = FALSE, 
    final.dens = FALSE, errfn = rsse, ...) 
{

    argl <- list(...)
    targind <- match(names(argl), "tol")
    targind <- which(!is.na(targind))
    margind <- match(names(argl), "method")
    margind <- which(!is.na(margind))
    if ((length(targind) == 0) & identical(abcmethod, abc)) {
        argl$tol <- 0.01
    }
    if ((length(margind) == 0) & identical(abcmethod, abc)) {
        argl$method <- "rejection"
    }

    if (!is.matrix(obs) | is.data.frame(obs)) {
        obs <- matrix(obs, nrow = 1)
    }
    if (!is.matrix(param) | is.data.frame(param)) {
        param <- as.matrix(param)
    }
    if (!is.matrix(sumstats) | is.data.frame(sumstats)) {
        sumstats <- as.matrix(sumstats)
    }
    if (!is.null(obspar) | is.data.frame(obspar)) {
        if (!is.matrix(obspar)) {
            obspar <- matrix(obspar, byrow = T, ncol = ncol(param))
        }
        if (nrow(obs) != nrow(obspar)) {
            stop("Please supply observed statistics and observed parameter matrices with the same number of rows!\n")
        }
    }
    argl$param=param

if (!length(colnames(param))) {
        colnames(param) <- paste("P", 1:ncol(param), sep = "")
}
if (!length(colnames(sumstats))) {
        colnames(sumstats) <- paste("C", 1:ncol(sumstats), sep = "")
}

    sumstats <- sumstats[, sumsubs]
    obs <- obs[, sumsubs]
    data <- matrix(obs, nrow = 1)
    npar <- ncol(param)
    nr <- nrow(sumstats)
    nstats <- length(sumsubs)
    q <- (!is.null(obspar)) & (do.err)
    if (!q) {
        do.err <- FALSE
    }
    cm <- combmat(nstats, limit)
    nc <- nrow(cm)

    dom<-(0%in%do.only)

    if(dom &!is.matrix(do.only)){
	do.only<-matrix(do.only,nrow=1)
    }

    if (is.null(do.only)|max(do.only)>nc) {
        do.only <- 1:nc
    }

    do.only<-as.matrix(do.only)

    if(dom & (ncol(do.only) > nstats)){
	stop("supplied subset matrix (do.only) does not match number of given summaries (ncol(sumsubs)) \n")
    }
    if(dom){
	ldo<-nrow(do.only)
    }
    else{
	ldo<-length(do.only)
    }

    vals <- err <- NULL
    critvals <- matrix(0, 1, length(do.only))
    best <- 1
#    for (i in 1:length(do.only)) {
    for (i in 1:ldo) {
	if(dom){
        	I <- which(do.only[i, ] == 1)
	}
	else{
        	I <- which(cm[do.only[i], ] == 1)
	}
        if (verbose) {
            cat("doing statistics: ", sumsubs[I], "   (", i, 
                "/", ldo, ")\n")
        }
	
	argl$target=data[I]
	argl$sumstat=sumstats[,I]
	valsi <- do.call(abcmethod,argl)
#     	valsi <- eval(parse(text=paste("abcmethod(data[I], param, sumstats[, I],",addstr,  "...)" )))

        if (is.null(valsi$adj.values)) {
            valsi <- valsi$unadj.values
        }
        else {
            valsi <- valsi$adj.values
        }
        if (do.err) {
            err[i] <- errfn(valsi, obspar, apply(param, 2, var))
        }
        if (do.crit) {
            cat("doing criterion calculation...\n")
            critvals[i] <- crit(valsi)
            stick <- (critvals[best] <= critvals[i])
            best <- ifelse(stick, best, i)
            if (final.dens) {
                if ((i == 1) | !stick) {
                  vals <- valsi
                }
            }
        }
    }

    err <- err[best]

    l <- list()
    if (do.crit) {
        l$critvals <- critvals
	if(dom){
		besti<-best
    		best <- matrix(which(do.only[besti, ]==1),nrow=1)
	}
	else{
        	besti <- do.only[best]
    		best <- matrix(which(cm[besti, ] == 1),nrow=1)
	}
	rownames(best)<-besti

        l$best <- best
    }
    if (do.err) {
        l$err <- err
    }
    if (final.dens) {
        l$post.sample <- vals
    }
    l$posssubs <- do.only
    l$sumsubs <- sumsubs
    return(l)
}
