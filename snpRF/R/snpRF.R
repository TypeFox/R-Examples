## mylevels() returns levels if given a factor, otherwise 0.
mylevels <- function(x) if (is.factor(x)) levels(x) else 0

snpRF <-
    function(x.autosome=NULL,x.xchrom=NULL, xchrom.names=NULL, x.covar=NULL, y,  
    	     xtest.autosome=NULL,xtest.xchrom=NULL, xtest.covar=NULL,
    	     ytest=NULL, ntree=500,
             mtry=floor(sqrt(sum(c(ncol(x.autosome),ncol(x.xchrom)/2,ncol(x.covar))))),
             replace=TRUE, classwt=NULL, cutoff, strata,
             sampsize = if (replace) max(c(nrow(x.autosome),nrow(x.xchrom),nrow(x.covar))) else ceiling(.632*max(c(nrow(x.autosome),nrow(x.xchrom),nrow(x.covar)))),
             nodesize = 1,
             maxnodes=NULL,
             importance=FALSE, localImp=FALSE,
             proximity, oob.prox=proximity,
             norm.votes=TRUE, do.trace=FALSE,
             keep.forest=!is.null(y) && (is.null(xtest.autosome) & is.null(xtest.xchrom) & is.null(xtest.covar)),
             keep.inbag=FALSE, ...) {
    addclass <- is.null(y)
    classRF <- addclass || is.factor(y)
    if (!classRF && length(unique(y)) <= 5) {
        warning("The response has five or fewer unique values.  Are you sure you want to do regression?")
    }
    if (classRF && !addclass && length(unique(y)) < 2)
        stop("Need at least two classes to do classification.")


    ### must have at least one of x.* ###

    if(is.null(x.autosome) & is.null(x.xchrom) & is.null(x.covar)) stop("must have one of x.autosome, x.xchrom, or x.covar")

    ### force predictor data to be matrices ###

    if(!is.null(x.autosome) & class(x.autosome)!="matrix") stop("x.autosome must be a matrix") 
    if(!is.null(x.xchrom) & class(x.xchrom)!="matrix") stop("x.xchrom must be a matrix") 
    if(!is.null(x.covar) & class(x.covar)!="matrix") stop("x.covar must be a matrix") 

    ### make sure that x.xchrom are all 0's and 1's ###

    if(!all(unique(as.vector(x.xchrom)) %in% c(0,1))) stop("x.xchrom must be a matrix of 0 & 1s")

    ### check to see whether the xchrom and autosome have same number of people ###

    n<-as.numeric(as.character(unique(c(nrow(x.autosome),nrow(x.xchrom),nrow(x.covar)))))

    if(length(n)>1) stop("different number of obs. for x objects")

    p.xchrom<-ifelse(is.null(x.xchrom),0,ncol(x.xchrom))
    p.autosome<-ifelse(is.null(x.autosome),0,ncol(x.autosome))
    p.covar<-ifelse(is.null(x.covar),0,ncol(x.covar))


    x.row.names <- list()
    
    x.row.names[[1]]<-rownames(x.autosome)
    x.row.names[[2]]<-rownames(x.xchrom)
    x.row.names[[3]]<-rownames(x.covar)

    x.row.names<-do.call(rbind,lapply(x.row.names,function(x) dimnames(x)[[2]]))

    x.row.names<-x.row.names[!duplicated(x.row.names),]

    if(!is.null(x.row.names) && dim(x.row.names)[1]>1) stop("row.names do not agree between x.*")
    
    if(!is.null(x.autosome) && is.null(dimnames(x.autosome)[[2]])) dimnames(x.autosome)[[2]]<-paste("A",seq(1,p.autosome),sep="") 

    if(!is.null(x.xchrom) && is.null(xchrom.names)) { xchrom.names<-paste("X",seq(1,p.xchrom/2),sep="")
    }else{ if(length(xchrom.names) != p.xchrom/2) stop("length(xchrom.names) != dim(x.xchrom)[2]/2")}

    if(!is.null(x.covar) && is.null(dimnames(x.covar)[[2]])) dimnames(x.covar)[[2]]<-paste("C",seq(1,p.covar),sep="") 

    x.col.names <- c(dimnames(x.autosome)[[2]],dimnames(x.covar)[[2]],xchrom.names)

    ## overcome R's lazy evaluation:
    keep.forest <- keep.forest

    testdat <- !is.null(xtest.autosome) | !is.null(xtest.xchrom) |  !is.null(xtest.covar)

    if (testdat) {
	if((!is.null(x.autosome) & is.null(xtest.autosome)) | 
	   (!is.null(x.xchrom) & is.null(xtest.xchrom)) | 
	   (!is.null(x.covar) & is.null(xtest.covar)))        stop("each x object must have a corresponding xtest object") 
    
        if ((!is.null(xtest.autosome) && ncol(x.autosome) != ncol(xtest.autosome)) | 
	    (!is.null(xtest.xchrom) && ncol(x.xchrom) != ncol(xtest.xchrom)) | 
	    (!is.null(xtest.covar) && ncol(x.covar) != ncol(xtest.covar))) stop("x and xtest must have same number of columns")

        ntest <- as.numeric(as.character(unique(c(nrow(xtest.autosome),nrow(xtest.xchrom),nrow(xtest.covar)))))
	if(length(ntest)>1) stop("xtest datasets must have same number of observations")

    	xts.row.names <- list()
    
	xts.row.names[[1]]<-rownames(xtest.autosome)
    	xts.row.names[[2]]<-rownames(xtest.xchrom)
    	xts.row.names[[3]]<-rownames(xtest.covar)

    	xts.row.names<-do.call(rbind,lapply(xts.row.names,function(x) dimnames(x)[[2]]))

    	xts.row.names<-xts.row.names[!duplicated(xts.row.names),]

    	if(!is.null(xts.row.names) && dim(xts.row.names)[1]>1) stop("row.names do not agree between xtest.*")
    
    }

    ## Make sure mtry is in reasonable range.
    if ((mtry < 1) || (mtry > (p.autosome + (p.xchrom/2) + p.covar)))
        warning("invalid mtry: reset to within valid range")
    mtry <- max(1, min((p.autosome + (p.xchrom/2) + p.covar), round(mtry)))


    ## Check for NAs.

    if ( (!is.null(x.autosome) && any(is.na(x.autosome))) | 
         (!is.null(x.xchrom) && any(is.na(x.xchrom))) | 
         (!is.null(x.covar) && any(is.na(x.covar)))) stop("NA not permitted in genetic predictors")

    if (testdat && 
       	((!is.null(xtest.autosome) && any(is.na(xtest.autosome))) | 
         (!is.null(xtest.xchrom) && any(is.na(xtest.xchrom))) |
	 (!is.null(xtest.covar) && any(is.na(xtest.covar)))))     stop("NA not permitted in genetic xtest")

    if (any(is.na(y))) stop("NA not permitted in response")
    if (!is.null(ytest) && any(is.na(ytest))) stop("NA not permitted in ytest")


    ### forcing predictor variables to be matrices ###


    ncat.xchrom <- rep(1, p.xchrom)
    ncat.autosome <- rep(1, p.autosome+p.covar)
    xlevels.xchrom <- as.list(rep(0, p.xchrom))
    xlevels.autosome <- as.list(rep(0, p.autosome))
    xlevels.covar <- as.list(rep(0, p.covar))

    maxcat<-1


    if (classRF) {
        nclass <- length(levels(y))
        ## Check for empty classes:
        if (any(table(y) == 0)) stop("Can't have empty classes in y.")
        if (!is.null(ytest)) {
            if (!is.factor(ytest)) stop("ytest must be a factor")
            if (!all(levels(y) == levels(ytest)))
                stop("y and ytest must have the same levels")
        }
        if (missing(cutoff)) {
            cutoff <- rep(1 / nclass, nclass)
        } else {
            if (sum(cutoff) > 1 || sum(cutoff) < 0 || !all(cutoff > 0) ||
                length(cutoff) != nclass) {
                stop("Incorrect cutoff specified.")
            }
            if (!is.null(names(cutoff))) {
                if (!all(names(cutoff) %in% levels(y))) {
                    stop("Wrong name(s) for cutoff")
                }
                cutoff <- cutoff[levels(y)]
            }
        }
        if (!is.null(classwt)) {
            if (length(classwt) != nclass)
                stop("length of classwt not equal to number of classes")
            ## If classwt has names, match to class labels.
            if (!is.null(names(classwt))) {
                if (!all(names(classwt) %in% levels(y))) {
                    stop("Wrong name(s) for classwt")
                }
                classwt <- classwt[levels(y)]
            }
            if (any(classwt <= 0)) stop("classwt must be positive")
            ipi <- 1
        } else {
            classwt <- rep(1, nclass)
            ipi <- 0
        }
    } else addclass <- FALSE

    if (missing(proximity)) proximity <- addclass
    if (proximity) {
        prox <- matrix(0.0, n, n)
        proxts <- if (testdat) matrix(0, ntest, ntest + n) else double(1)
    } else {
        prox <- proxts <- double(1)
    }

    if (localImp) {
        importance <- TRUE
        impmat <- matrix(0, ((p.xchrom/2)+p.autosome+p.covar), n)
    } else impmat <- double(1)

    if (importance) {
        #if (nPerm < 1) nPerm <- as.integer(1) else nPerm <- as.integer(nPerm)
        if (classRF) {

            impout <- matrix(0.0, ((p.xchrom/2)+p.autosome+p.covar), nclass + 2)
            impSD <- matrix(0.0,((p.xchrom/2)+p.autosome+p.covar) , nclass + 1)
        } else {
            impout <- matrix(0.0, ((p.xchrom/2)+p.autosome+p.covar), 2)
            impSD <- double(((p.xchrom/2)+p.autosome+p.covar))
            names(impSD) <- x.col.names
        }
    } else {
        impout <- double(((p.xchrom/2)+p.autosome+p.covar))
        impSD <- double(1)
    }

    nsample <- if (addclass) 2 * n else n
    Stratify <- length(sampsize) > 1
    if ((!Stratify) && sampsize > n) stop("sampsize too large")
    if (Stratify && (!classRF)) stop("sampsize should be of length one")
    if (classRF) {
        if (Stratify) {
            if (missing(strata)) strata <- y
            if (!is.factor(strata)) strata <- as.factor(strata)
            nsum <- sum(sampsize)
            if (length(sampsize) > nlevels(strata))
                stop("sampsize has too many elements.")
            if (any(sampsize <= 0) || nsum == 0)
                stop("Bad sampsize specification")
            ## If sampsize has names, match to class labels.
            if (!is.null(names(sampsize))) {
                sampsize <- sampsize[levels(strata)]
            }
            if (any(sampsize > table(strata)))
              stop("sampsize can not be larger than class frequency")
        } else {
            nsum <- sampsize
        }
        nrnodes <- 2 * trunc(nsum / nodesize) + 1
    } else {
        ## For regression trees, need to do this to get maximal trees.
        nrnodes <- 2 * trunc(sampsize/max(1, nodesize - 4)) + 1
    }
    if (!is.null(maxnodes)) {
        ## convert # of terminal nodes to total # of nodes
        maxnodes <- 2 * maxnodes - 1
        if (maxnodes > nrnodes) warning("maxnodes exceeds its max value.")
        nrnodes <- min(c(nrnodes, max(c(maxnodes, 1))))
    }


    #### combine x.autosome and x.covar for simplicity, since no special care ###
    #### needs to be taken for either of these sets of variables ###

    x.autosome<-cbind(x.autosome,x.covar)

    ## Compiled code expects variables in rows and observations in columns.

    ## in case x.* is NULL set as 0, this will work because x.* will have dim 0 ##
    ## and will be skipped in rf.c ##

    if(is.null(x.autosome)) {
      x.autosome<-0
    }else{
      x.autosome <- t(x.autosome)
    }
    if(is.null(x.xchrom)) {
      x.xchrom<-0
    }else{
      x.xchrom <- t(x.xchrom)
    }
    
    storage.mode(x.autosome) <- "double"
    storage.mode(x.xchrom) <- "double"

    if (testdat) {

    	#### combine xtest.autosome and xtest.covar for simplicity, as above for ###
    	#### x.autosome and x.covar ###

	xtest.autosome<-cbind(xtest.autosome,xtest.covar)

        if(is.null(xtest.autosome)) {
	  xtest.autosome<-0
	  p.xta<-0
	}else{
	  p.xta<-dim(xtest.autosome)[2]
	  xtest.autosome <- t(xtest.autosome)
        }        
	if(is.null(xtest.xchrom)) {
	  xtest.xchrom<-0
	  p.xtx<-0

	}else{
	  p.xtx<-dim(xtest.xchrom)[2]
	  xtest.xchrom <- t(xtest.xchrom)
        }
        storage.mode(xtest.autosome) <- "double"
        storage.mode(xtest.xchrom) <- "double"

        if (is.null(ytest)) {
            ytest <- labelts <- 0
        } else {
            labelts <- TRUE
        }
    } else {
        xtest.autosome <- double(1)
	p.xta<-0
        xtest.xchrom <- double(1)
	p.xtx<-0
        ytest <- double(1)
        ntest <- 1
        labelts <- FALSE
    }
    nt <- if (keep.forest) ntree else 1

    if (classRF) {
        cwt <- classwt
        threshold <- cutoff
        error.test <- if (labelts) double((nclass+1) * ntree) else double(1)

        rfout <- .C("classRF",

                    autosome=x.autosome,
		    xchrom=x.xchrom,
                    xdim = as.integer(c(p.autosome+p.covar,p.xchrom,n)),

                    y = as.integer(y),
                    nclass = as.integer(nclass),

                    ncat = as.integer(c(ncat.autosome,ncat.xchrom)),

                    maxcat = as.integer(maxcat),
                    sampsize = as.integer(sampsize),
                    strata = if (Stratify) as.integer(strata) else integer(1),
                    Options = as.integer(c(addclass,
                    importance,
                    localImp,
                    proximity,
                    oob.prox,
                    do.trace,
                    keep.forest,
                    replace,
                    Stratify,
                    keep.inbag)),
                    ntree = as.integer(ntree),
                    mtry = as.integer(mtry),
                    ipi = as.integer(ipi),
                    classwt = as.double(cwt),
                    cutoff = as.double(threshold),
                    nodesize = as.integer(nodesize),
                    outcl = integer(nsample),
                    counttr = integer(nclass * nsample),
                    prox = prox,
                    impout = impout,
                    impSD = impSD,
                    impmat = impmat,
                    nrnodes = as.integer(nrnodes),
                    ndbigtree = integer(ntree),
                    nodestatus = integer(nt * nrnodes),
                    bestvar = integer(nt * nrnodes),
                    treemap = integer(nt * 2 * nrnodes),
                    nodepred = integer(nt * nrnodes),
                    xbestsplit = double(nt * nrnodes),
                    errtr = double((nclass+1) * ntree),
                    testdat = as.integer(testdat),

                    xtsA = as.double(xtest.autosome),
		    xtsX = as.double(xtest.xchrom),
		    xtsDim = as.integer(c(p.xta,p.xtx)),

                    clts = as.integer(ytest),
                    nts = as.integer(ntest),
                    countts = double(nclass * ntest),
                    outclts = as.integer(numeric(ntest)),
                    labelts = as.integer(labelts),
                    proxts = proxts,
                    errts = error.test,
                    inbag = if (keep.inbag)
                    matrix(integer(n * ntree), n) else integer(n),
                    PACKAGE="snpRF")[-1]


        if (keep.forest) {
            ## deal with the random forest outputs
            max.nodes <- max(rfout$ndbigtree)
            treemap <- aperm(array(rfout$treemap, dim = c(2, nrnodes, ntree)),
                             c(2, 1, 3))[1:max.nodes, , , drop=FALSE]
        }
        if (!addclass) {
            ## Turn the predicted class into a factor like y.
            out.class <- factor(rfout$outcl, levels=1:nclass,
                                labels=levels(y))
            names(out.class) <- x.row.names
            con <- table(observed = y,
                         predicted = out.class)[levels(y), levels(y)]
            con <- cbind(con, class.error = 1 - diag(con)/rowSums(con))
        }
        out.votes <- t(matrix(rfout$counttr, nclass, nsample))[1:n, ]
        oob.times <- rowSums(out.votes)
        if (norm.votes)
            out.votes <- t(apply(out.votes, 1, function(x) x/sum(x)))

        dimnames(out.votes) <- list(x.row.names, levels(y))

        class(out.votes) <- c(class(out.votes), "votes")
        if (testdat) {
            out.class.ts <- factor(rfout$outclts, levels=1:nclass,
                                   labels=levels(y))
            names(out.class.ts) <- xts.row.names
            out.votes.ts <- t(matrix(rfout$countts, nclass, ntest))
            dimnames(out.votes.ts) <- list(xts.row.names, levels(y))
            if (norm.votes)
                out.votes.ts <- t(apply(out.votes.ts, 1,
                                        function(x) x/sum(x)))
            class(out.votes.ts) <- c(class(out.votes.ts), "votes")
            if (labelts) {
                testcon <- table(observed = ytest,
                                 predicted = out.class.ts)[levels(y), levels(y)]
                testcon <- cbind(testcon,
                                 class.error = 1 - diag(testcon)/rowSums(testcon))
            }
        }
        cl <- match.call()
        cl[[1]] <- as.name("snpRF")

        out <- list(call = cl,
                    type = if (addclass) "unsupervised" else "classification",
                    predicted = if (addclass) NULL else out.class,
                    err.rate = if (addclass) NULL else t(matrix(rfout$errtr,
                    nclass+1, ntree,
                    dimnames=list(c("OOB", levels(y)), NULL))),
                    confusion = if (addclass) NULL else con,
                    votes = out.votes,
                    oob.times = oob.times,
                    classes = levels(y),
                    importance = if (importance)
                    matrix(rfout$impout, (p.xchrom/2) + p.autosome + p.covar, nclass+2,
                           dimnames = list(x.col.names,
                           c(levels(y), "MeanDecreaseAccuracy",
                             "MeanDecreaseGini")))
                    else matrix(rfout$impout, ncol=1,
                                dimnames=list(x.col.names, "MeanDecreaseGini")),
                    importanceSD = if (importance)
                    matrix(rfout$impSD, (p.xchrom/2)+p.autosome+p.covar, nclass + 1,
                           dimnames = list(x.col.names,
                           c(levels(y), "MeanDecreaseAccuracy")))
                    else NULL,
                    localImportance = if (localImp)
                    matrix(rfout$impmat, (p.xchrom/2)+p.autosome+p.covar, n,
                           dimnames = list(x.col.names,x.row.names)) else NULL,
                    proximity = if (proximity) matrix(rfout$prox, n, n,
                    dimnames = list(x.row.names, x.row.names)) else NULL,
                    ntree = ntree,
                    mtry = mtry,
                    forest = if (!keep.forest) NULL else {
                        list(ndbigtree = rfout$ndbigtree,
                             nodestatus = matrix(rfout$nodestatus,
                             ncol = ntree)[1:max.nodes,, drop=FALSE],
                             bestvar = matrix(rfout$bestvar, ncol = ntree)[1:max.nodes,, drop=FALSE],  ### this is after X markers are converted to 1 column
                             treemap = treemap,
                             nodepred = matrix(rfout$nodepred,
                             ncol = ntree)[1:max.nodes,, drop=FALSE],
                             xbestsplit = matrix(rfout$xbestsplit,
                             ncol = ntree)[1:max.nodes,, drop=FALSE],
                             pid = rfout$classwt, cutoff=cutoff, ncat=c(ncat.autosome,ncat.xchrom),  ### this counts X markers twice
			     ncat.autosome=ncat.autosome,ncat.xchrom=ncat.xchrom, ### this counts X markers twice, 
			     							  ### and splits up autosomal+covar and X markers
                             maxcat = maxcat,
                             nrnodes = max.nodes, ntree = ntree,
                             nclass = nclass, xlevels=c(xlevels.autosome,xlevels.covar,xlevels.xchrom))
                    },
                    y = if (addclass) NULL else y,
                    test = if(!testdat) NULL else list(
                    predicted = out.class.ts,
                    err.rate = if (labelts) t(matrix(rfout$errts, nclass+1,
                    ntree,
                dimnames=list(c("Test", levels(y)), NULL))) else NULL,
                    confusion = if (labelts) testcon else NULL,
                    votes = out.votes.ts,
                    proximity = if(proximity) matrix(rfout$proxts, nrow=ntest,
                    dimnames = list(xts.row.names, c(xts.row.names,
                    x.row.names))) else NULL),
                    inbag = if (keep.inbag) rfout$inbag else NULL)
    } else {
    
	stop("Regression trees not implemented")
	out<-NULL

    }
    class(out) <- "snpRF"
    return(out)
}
