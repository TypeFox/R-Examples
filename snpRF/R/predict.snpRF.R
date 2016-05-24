"predict.snpRF" <-
    function (object, newdata.autosome=NULL, newdata.xchrom=NULL, xchrom.names=NULL, 
    	      newdata.covar=NULL, type = "response", norm.votes = TRUE,
              predict.all=FALSE, proximity = FALSE, nodes=FALSE, cutoff, ...)
{
    if (!inherits(object, "snpRF"))
        stop("object not of class snpRF")
    if (is.null(object$forest)) stop("No forest component in the object")
    out.type <- charmatch(tolower(type),
                          c("response", "prob", "vote", "class"))
    if (is.na(out.type))
        stop("type must be one of 'response', 'prob', 'vote'")
    if (out.type == 4) out.type <- 1
    if (out.type != 1 && object$type == "regression")
        stop("'prob' or 'vote' not meaningful for regression")
    if (out.type == 2)
        norm.votes <- TRUE
    if (is.null(newdata.autosome) & is.null(newdata.autosome) & is.null(newdata.autosome)) {
        if (object$type == "regression") return(object$predicted)
        if (proximity & is.null(object$proximity))
            warning("cannot return proximity without new data if random forest object does not already have proximity")
        if (out.type == 1) {
            if (proximity) {
                return(list(pred = object$predicted,
                            proximity = object$proximity))
            } else return(object$predicted)
        }
        if (norm.votes) {
            t1 <- t(apply(object$votes, 1, function(x) { x/sum(x) }))
            class(t1) <- c(class(t1), "votes")
            if (proximity) return(list(pred = t1, proximity = object$proximity))
            else return(t1)
        } else {
            if (proximity) return(list(pred = object$votes, proximity = object$proximity))
            else return(object$votes)
        }
    }
    if (missing(cutoff)) {
        cutoff <- object$forest$cutoff
    } else {
        if (sum(cutoff) > 1 || sum(cutoff) < 0 || !all(cutoff > 0) ||
            length(cutoff) != length(object$classes)) {
            stop("Incorrect cutoff specified.")
        }
        if (!is.null(names(cutoff))) {
            if (!all(names(cutoff) %in% object$classes)) {
                stop("Wrong name(s) for cutoff")
            }
            cutoff <- cutoff[object$classes]
        }
    }

    if (object$type == "unsupervised")
        stop("Can't predict unsupervised forest.")

    #if (inherits(object, "randomForest.formula")) {
    #    newdata <- as.data.frame(newdata)
    #    rn <- row.names(newdata)
    #    Terms <- delete.response(object$terms)
    #    x <- model.frame(Terms, newdata, na.action = na.omit)
    #    keep <- match(row.names(x), rn)
    #} else {
        
	x <- list(newdata.autosome,newdata.covar,newdata.xchrom)
	rm(newdata.autosome,newdata.covar,newdata.xchrom)

	class.x<-unlist(lapply(x,class))
	if(any(class.x[class.x!="NULL"]!="matrix")) stop("newdata.* must be matrices")

	n.p<-unlist(lapply(x,function(z) ifelse(is.null(z),0,nrow(z))))

	if(length(unique(n.p[n.p!=0]))!=1) stop("newdata.* have different no. rows")

	n.p<-as.numeric(as.character(unique(n.p[n.p!=0])))

        if (any(unlist(lapply(x[class.x!=NULL],is.na)))) 
            stop("missing values in newdata.*")
        keep <- 1:n.p
	     
	rn<-do.call(rbind,lapply(x,function(z) dimnames(z)[[1]]))
	rn<-rn[!duplicated(rn),,drop=F]

	if(!is.null(rn) && dim(rn)[1]>1) stop("differing row names for newdata.*")
	
        if (is.null(rn)) { rn <- keep }else{rn<-as.vector(rn)}

    #}
    vname <- if (is.null(dim(object$importance))) {
        names(object$importance)
    } else {
        rownames(object$importance)
    }

    if(!is.null(x[[1]]) && is.null(dimnames(x[[1]])[[2]])) dimnames(x[[1]])[[2]]<-paste("A",1:dim(x[[1]])[2],sep="")
    if(!is.null(x[[2]]) && is.null(dimnames(x[[2]])[[2]])) dimnames(x[[2]])[[2]]<-paste("C",1:dim(x[[2]])[2],sep="")
    if(!is.null(x[[3]]) && is.null(xchrom.names)) xchrom.names<-paste("X",1:(dim(x[[3]])[2]/2),sep="")

    newdata.vname<-list(dimnames(x[[1]])[[2]],dimnames(x[[2]])[[2]],xchrom.names)
    
    if (length(unlist(newdata.vname)) != length(vname)) {

       stop("number of variables in newdata.* does not match that in the training data")
    
    } else {
        if (any(! (vname %in% unlist(newdata.vname)))) stop("variables in the training data missing in newdata")

	n.c<-cumsum(unlist(lapply(x,function(z) ifelse(is.null(z),0,dim(z)[2]))))

	if(!is.null(x[[1]])) x[[1]] <- x[[1]][, vname[vname %in% newdata.vname[[1]]], drop=FALSE]

        if(!is.null(x[[2]])) x[[2]] <- x[[2]][, vname[vname %in% newdata.vname[[2]]], drop=FALSE]

	if(!is.null(x[[3]])) { 

			     tmp<-as.numeric(rep(NA,length(newdata.vname[[3]])*2))
			     tmp[seq(1,length(tmp)-1,2)]<-(2*match(newdata.vname[[3]],vname))-1
			     tmp[seq(2,length(tmp),2)]<-2*match(newdata.vname[[3]],vname)
        		     x[[3]] <- x[[3]][, tmp-(2*n.c[2]), drop=FALSE]
			     rm(tmp)
	}

    }

    #### combine x.autosome and x.covar for simplicity, since no special care ###
    #### needs to be taken for either of these sets of variables ###

    x.autosome<-cbind(x[[1]],x[[2]])
    x.xchrom<-x[[3]]
    rm(x)

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
    

    ### data.frames not allowed ###
    #if (is.data.frame(x)) {
    #    xfactor <- which(sapply(x, is.factor))
    #    if (length(xfactor) > 0 && "xlevels" %in% names(object$forest)) {
    #        for (i in xfactor) {
    #            if (any(! levels(x[[i]]) %in% object$forest$xlevels[[i]]))
    #                stop("New factor levels not present in the training data")
    #            x[[i]] <-
    #                factor(x[[i]],
    #                       levels=levels(x[[i]])[match(levels(x[[i]]), object$forest$xlevels[[i]])])
    #        }
    #    }
    #    cat.new <- sapply(x, function(x) if (is.factor(x) && !is.ordered(x))
    #                      length(levels(x)) else 1)
    #    if (!all(object$forest$ncat == cat.new))
    #        stop("Type of predictors in new data do not match that of the training data.")
    #}

    mdim <- c(nrow(x.autosome),nrow(x.xchrom))
    ntest <- n.p
    ntree <- object$forest$ntree
    maxcat <- max(object$forest$ncat)
    nclass <- object$forest$nclass
    nrnodes <- object$forest$nrnodes
    ## get rid of warning:
    op <- options(warn=-1)
    on.exit(options(op))

    if (predict.all) {
       	    treepred <- matrix(integer(ntest * ntree), ncol=ntree)
    }else{
        treepred <- numeric(ntest)
    }
    proxmatrix <- if (proximity) matrix(0, ntest, ntest) else numeric(1)
    nodexts <- if (nodes) integer(ntest * ntree) else integer(ntest)

    countts <- matrix(0, ntest, nclass)


    t1 <- .C("classForest",
             xdim = as.integer(c(mdim,ntest)),
	     autosome=x.autosome,
	     xchrom=x.xchrom,
             nclass = as.integer(object$forest$nclass),
             maxcat = as.integer(maxcat),
             nrnodes = as.integer(nrnodes),
             jbt = as.integer(ntree),
             xbestsplit = as.double(object$forest$xbestsplit),
             pid = object$forest$pid,
             cutoff = as.double(cutoff),
             countts = as.double(countts),
             treemap = as.integer(aperm(object$forest$treemap,
                                 c(2, 1, 3))),
             nodestatus = as.integer(object$forest$nodestatus),
             cat = as.integer(object$forest$ncat),
             nodepred = as.integer(object$forest$nodepred),
             treepred = as.integer(treepred),
             jet = as.integer(numeric(ntest)),
             bestvar = as.integer(object$forest$bestvar),
             nodexts = as.integer(nodexts),
             ndbigtree = as.integer(object$forest$ndbigtree),
             predict.all = as.integer(predict.all),
             prox = as.integer(proximity),
             proxmatrix = as.double(proxmatrix),
             nodes = as.integer(nodes),
             PACKAGE = "snpRF")


        if (out.type > 1) {
            out.class.votes <- t(matrix(t1$countts, nrow = nclass, ncol = ntest))
            if (norm.votes)
                out.class.votes <-
                    sweep(out.class.votes, 1, rowSums(out.class.votes), "/")
            z <- matrix(NA, length(rn), nclass,
                        dimnames=list(rn, object$classes))
            z[keep, ] <- out.class.votes
             class(z) <- c(class(z), "votes")
            res <- z
        } else {
            out.class <- factor(rep(NA, length(rn)),
                                levels=1:length(object$classes),
                                labels=object$classes)
            out.class[keep] <- object$classes[t1$jet]
            names(out.class)[keep] <- rn[keep]
            res <- out.class
        }
        if (predict.all) {
            treepred <- matrix(object$classes[t1$treepred],
                               nrow=length(keep), dimnames=list(rn[keep], NULL))
            res <- list(aggregate=res, individual=treepred)
        }
        if (proximity)
            res <- list(predicted = res, proximity = structure(t1$proxmatrix,
                                         dim = c(ntest, ntest),
                                         dimnames = list(rn[keep], rn[keep])))
        if (nodes) attr(res, "nodes") <- matrix(t1$nodexts, ntest, ntree,
                                                dimnames=list(rn[keep], 1:ntree))
    
    res
}
