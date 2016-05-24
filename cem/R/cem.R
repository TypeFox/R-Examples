group.var <- function(x, groups){ 
	n <- length(x)
#	print(str(x))
	tmp <- numeric(n)
	ngr <- length(groups)
	all.idx <- NULL
#	cat(sprintf("n=%d, ngr=%d\n", n, ngr))
	 for(i in 1:ngr){
	  idx <- which(x %in% groups[[i]])
	   if(length(idx)>0){
		tmp[idx] <- sprintf("g%.2d",i) 
		all.idx <- c(all.idx, idx)
	   }
	  } 
	  if(length(all.idx)>0 && length(all.idx)<n)
	   tmp[-all.idx] <- x[-all.idx]
	tmp <- factor(tmp)
}

`cem` <-
function (treatment=NULL, data = NULL, datalist=NULL, cutpoints = NULL,  
    grouping = NULL, drop=NULL, eval.imbalance = FALSE, k2k=FALSE,  
	method=NULL, mpower=2, L1.breaks = NULL, L1.grouping = NULL, 
    verbose = 0, baseline.group="1",keep.all=FALSE)
{
    L1data <- data
	L1datalist <- datalist
    if(k2k)  keep.all <- TRUE
	if(!is.null(grouping) & !is.null(names(grouping))){
      gn <- names(grouping)
      n.gn <- length(gn)
      nd <- 0
      if(!is.null(datalist) & is.list(datalist))
       nd <- length(datalist)
	  for(g in 1:n.gn){
	   if(nd>0){
        for(d in 1:nd)
		 datalist[[d]][[gn[g]]]  <- group.var( datalist[[d]][[gn[g]]] , grouping[[g]])
	   }
	   if(!is.null(data))
	    data[[gn[g]]] <- group.var(data[[gn[g]]], grouping[[g]])
	   if(!is.null(cutpoints[gn[g]]))
	    cutpoints[gn[g]] <- NULL   			   
	  }
	 }
   if(!is.null(data) & is.null(datalist)){
     mat <- cem.main(treatment=treatment, data=data, cutpoints = cutpoints,  drop=drop, 
               k2k=k2k, method=method, mpower=mpower, 
			   verbose = verbose, baseline.group=baseline.group,
       keep.all=keep.all)
 	 if(eval.imbalance & !is.null(treatment) & length(which(mat$matched==TRUE)>0))
	  mat$imbalance <- imbalance(mat$groups, data=L1data, drop=mat$drop, breaks=L1.breaks,
	                  weights=mat$w, grouping=L1.grouping)
	mat$grouping <- grouping				  
    return( mat )
   } 


   if (is.null(datalist) & is.null(data)) 
	 stop("Specify at list `data' argument", call. = FALSE)
    if (is.data.frame(datalist)) {
        stop("Please specify a list of data frames", call. = FALSE)
    }
    if (!is.null(datalist) & !is.list(datalist)) {
        stop("Please specify a list of data frames", call. = FALSE)
    }

   nd <- length(datalist)
   list.obj <- vector(nd, mode="list")


   if(!is.null(data)){
	n <- dim(data)[1]
    multi.strata <- matrix(NA, n, nd)
        	
    tmp <- na.omit(data)
    comp <- match(rownames(tmp), rownames(data))
	mis <- 1:dim(data)[1] 
	mis <- mis[-comp]
	n.mis <- length(mis)
    n.comp <- length(comp)
	bigdata <- NULL
	
	for(i in 1:nd)
	   bigdata <- rbind(bigdata, datalist[[i]])

    rownames(bigdata) <- 1:(dim(bigdata)[1])
	mat <- cem.main(treatment=treatment, data=bigdata, cutpoints = cutpoints,  drop=drop, 
               k2k=k2k, method=method, mpower=mpower, 
			   verbose = verbose, baseline.group=baseline.group,keep.all=keep.all)
	
	for(i in 1:nd){
     all <- c(1:n + (i-1)*n )
	 multi.strata[,i] <- mat$strata[all]
    }
	g <- function(x){
	 tmp <- table(x)
	 id <- as.integer(which.max(tmp))
	 names(tmp)[id]
	}
	
	new.strata <- as.integer(apply(multi.strata, 1, g))

    for(i in 1){
     obj <- mat
	 all <- c(1:n + (i-1)*n )
    if(keep.all){ 
        obj$X <- mat$X[all,]
	rownames(obj$X) <- rownames(datalist[[i]]) 
       }
    obj$strata <- new.strata
    obj$drop <- mat$drop
    obj$breaks <- mat$breaks
    imbalance <- NULL
    tab <- NULL
    obj$treatment <- treatment
    obj$n <- n
    if (!is.null(treatment)) {
        obj$groups <- mat$groups[all]
        obj$g.names <- levels(obj$groups)
        obj$n.groups <- length(obj$g.names)
        obj$group.idx <- sapply(obj$g.names, function(x) which(obj$groups == 
            x))
        names(obj$group.idx) <- paste("G", obj$g.names, sep = "")
        obj$group.len <- unlist(lapply(obj$group.idx, length))
        tmp <- find.strata(obj)
        obj$mstrata <- tmp$mstrata
        obj$mstrataID <- tmp$mstrataID
        obj$matched <- !is.na(obj$mstrata)
    }
    obj$imbalance <- imbalance
    obj$k2k <- k2k
    obj$w <- cem.weights(obj)
    class(obj) <- "cem.match"
    if (k2k) 
        obj <- k2k(obj, data=datalist[[i]], method = method, mpower = mpower, 
            verbose = verbose)
    if (!is.null(treatment)) 
        tab <- cem.summary(obj = obj, verbose = verbose)
    obj$tab <- tab
  	 
	 list.obj[[i]] <- obj
    }


    for(i in 2:nd){
     obj <- mat
	 all <- c(1:n + (i-1)*n )
    if(keep.all){
        obj$X <- mat$X[all,]
	rownames(obj$X) <- rownames(datalist[[i]])
    }
    obj$strata <- new.strata
    obj$drop <- mat$drop
    obj$breaks <- mat$breaks
    imbalance <- NULL
    tab <- NULL
    obj$treatment <- treatment
    obj$n <- n
    if (!is.null(treatment)) {
        obj$groups <- mat$groups[all]
        obj$g.names <- levels(obj$groups)
        obj$n.groups <- length(obj$g.names)
        obj$group.idx <- sapply(obj$g.names, function(x) which(obj$groups == 
            x))
        names(obj$group.idx) <- paste("G", obj$g.names, sep = "")
        obj$group.len <- unlist(lapply(obj$group.idx, length))
        tmp <- list.obj[[1]]
        obj$mstrata <- tmp$mstrata
        obj$mstrataID <- tmp$mstrataID
        obj$matched <- tmp$matched
    }
    obj$imbalance <- imbalance
    obj$k2k <- k2k
    obj$w <- cem.weights(obj)
    class(obj) <- "cem.match"
#    if (k2k) 
#        obj <- k2k(obj, data=datalist[[i]], method = method, mpower = mpower, 
#            verbose = verbose)
    if (!is.null(treatment)) 
        tab <- cem.summary(obj = obj, verbose = verbose)
    obj$tab <- tab
  	 
	 list.obj[[i]] <- obj
    }


	if(eval.imbalance  & !is.null(treatment)){
    avg.data <- datalist[[1]]
    nv <- length(colnames(datalist[[1]]))
    for(j in 1:nv){ 
	 if(is.numeric(avg.data[,j]) || is.integer(avg.data[,j])){
	  for(i in 2:nd)
       avg.data[,j] <- avg.data[,j] + L1datalist[[i]][,j]
	  avg.data[,j] <- avg.data[,j]/nd
	 } else {
	  tmp <- NULL
      for(i in 2:nd)
	   tmp <- cbind(tmp, L1datalist[[i]][,j])
	   avg.data[,j] <- factor(apply(tmp, 1, function(x) x[which.max(table(x))[1]])) 
	 }
     }
	 
    tmp <- imbalance(	list.obj[[1]]$groups, data=avg.data, drop=list.obj[[1]]$drop, breaks=L1.breaks,
	                  weights=list.obj[[1]]$w, grouping=L1.grouping) 
	for(i in 1:nd){
     list.obj[[i]]$imbalance <- tmp
	 list.obj[[i]]$grouping <- grouping
	}
    }
    unique <- TRUE
	
   } else {	
   for(i in 1:nd){
    list.obj[[i]] <- cem.main(treatment=treatment, data=datalist[[i]], cutpoints = cutpoints,  drop=drop, 
            k2k=k2k, method=method, mpower=mpower, verbose = verbose,
              baseline.group=baseline.group,keep.all=keep.all)
    unique <- FALSE
	if(eval.imbalance & !is.null(treatment))
	  list.obj[[i]]$imbalance <- imbalance(list.obj[[i]]$groups, data=L1datalist[[i]], drop=list.obj[[i]]$drop, breaks=L1.breaks,
	                  weights=list.obj[[i]]$w, grouping=L1.grouping)
    list.obj[[i]]$grouping <- grouping

    }
  
  }
  
   names(list.obj) <- sprintf("match%d",1:nd)
   class(list.obj) <- c("cem.match.list", "list")
   list.obj$unique <- unique
    
  return(list.obj)
}


print.cem.match <- function(x,...){
 if(!is.null(x$tab)){
  print(x$tab)
  cat("\n")
  }
 if(!is.null(x$imbalance))
  print(x$imbalance)
}


print.cem.match.list <- function(x,...){
 if(x$unique){
  if(!is.null(x[[1]]$tab)){
   print(x[[1]]$tab)
   cat("\n")
  }
  if(!is.null(x[[1]]$imbalance))
   print(x[[1]]$imbalance)
 } else {
  for(i in 1:(length(x)-1))
   print(x[i])
 }
}

