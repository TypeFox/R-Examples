`cem.main` <-
function (treatment=NULL, data, cutpoints = NULL,  drop=NULL, 
    k2k=FALSE, method=NULL, mpower=2, verbose = 0,
    baseline.group,keep.all=TRUE)
{
   if (is.null(data)) 
        stop("Dataframe must be specified", call. = FALSE)
    if (!is.data.frame(data)) {
        stop("Data must be a dataframe", call. = FALSE)
    }
      
	if(!is.null(treatment)){
	 groups <- as.factor(data[[treatment]])
	 drop <- c(drop,treatment)
	}
	drop <- unique(drop)
	dropped <- match(drop, colnames(data))
	dropped <- dropped[!is.na(dropped)]

	if(length(dropped)>0) 
		data <- data[-dropped]
	vnames <- colnames(data)
    if (sum(is.na(data)) > 0) 
        cat("The data contain missing values. CEM will match on them; see the manual for other options.\n")

    n <- dim(data)[1]
    nv <- dim(data)[2]
	mycut <- vector(nv, mode="list")
    names(mycut) <- vnames
	# preprocessing
	if(verbose > 1)
 	 cat("\npre-processing data")
	for (i in 1:nv) {	
	    if(verbose>1)
		 cat(".")
		tmp <- reduce.var(data[[i]], cutpoints[[vnames[i]]])
		data[[i]] <- tmp$x
		mycut[[vnames[i]]] <- tmp$breaks
    }

	obj <- cem.match(data = data, verbose = verbose)
	if(keep.all)
     obj$X <- data
	obj$drop <- drop
	obj$breaks <- mycut
	imbalance <- NULL
	tab <- NULL
	
	obj$treatment <- treatment
	obj$n <- dim(data)[1]

	if(!is.null(treatment)){
		obj$groups <- groups
		obj$g.names <- levels(obj$groups)
		obj$n.groups <- length(obj$g.names)
		obj$group.idx <- sapply(obj$g.names, function(x) which(obj$groups==x))
		names(obj$group.idx) <- paste("G",obj$g.names,sep="")
		obj$group.len <- unlist(lapply(obj$group.idx, length))
        
		tmp <- find.strata(obj)
		
		obj$mstrata <- tmp$mstrata
		obj$mstrataID <- tmp$mstrataID
		obj$matched <- !is.na(obj$mstrata)
        if(missing(baseline.group))
         baseline.group = "1"
        if(!(baseline.group %in% obj$g.names)| is.null(baseline.group)){ 
          if("1" %in% obj$g.names)
            baseline.group <- "1"
          else
            baseline.group <- obj$g.names[1]
        }
        if(verbose>0)
         cat(sprintf("\nUsing '%s'='%s' as baseline group\n", treatment, baseline.group)) 
        obj$baseline.group <- baseline.group		
	}

	if(!is.null(treatment))
	 tab <- cem.summary(obj=obj, verbose = verbose)

  
	obj$tab <- tab
    obj$k2k <- k2k
	obj$w <- cem.weights(obj)
	class(obj) <- "cem.match"
	if(k2k)
	 obj <- k2k(obj, data, method=method, mpower=mpower, verbose=verbose)
	return(invisible(obj))
}


cem.weights <- function (obj) 
{
    bg <- which(obj$g.names==obj$baseline.group)
    bgn <- sprintf("G%s",obj$baseline.group)
    
    w <- rep(0, obj$n)
    if (!is.null(obj$treatment)) {
        tmp <- table(obj$mstrata, obj$groups)
        wh <- t((sapply(1:NROW(tmp), function(x) tmp[x,bg]/tmp[x,]) * (obj$tab["Matched",]/obj$tab["Matched", bgn]))) 
        rownames(wh) <- rownames(tmp)
        colnames(wh) <- colnames(tmp)
        w <- numeric(obj$n)
        mID <- as.numeric(rownames(wh))
        nID <- length(mID)
        for(i in 1:nID){
            for(j in levels(obj$groups)){
                idx <- which(obj$mstrata==mID[i] & obj$groups==j)
                w[idx] <- wh[i, match(j, colnames(wh))]
            }
        }
    }
    w
}


old.cem.weights <- function(obj){
    w <- rep(0,obj$n)
	if(!is.null(obj$treatment)){
	 tmp <- table(obj$mstrata, obj$groups)
	 wh <- tmp[,2]/tmp[,1] * obj$tab[2,1]/obj$tab[2,2]
	 mID <- as.numeric(names(wh))
     idx <- match(obj$mstrata, mID)
	 idx2 <- match(mID[idx], mID)
     w <- as.numeric(wh[idx2])
	 w[which(is.na(w))] <- 0
	 w[which(obj$matched & (obj$groups ==  obj$g.names[2]))] <- 1	
	}
	w
}


find.strata <- function(obj){ 
 y <- unique(obj$strata)
 y <- y[!is.na(y)]
 mstrata <- match(obj$strata,y)
 n.st <- length(y)
 tab <- table(mstrata, obj$groups)
 tt <- apply(tab,1, function(x) all(x>0))
 idx <- which(tt == FALSE)
#idx <- as.integer(which(tab[,1]*tab[,2]<1)) 
 idx <- which(mstrata %in% idx)
 mstrata[idx] <- NA
 list(mstrata=mstrata, mstrataID=unique(na.omit(mstrata)))
}

