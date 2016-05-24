L1.meas <- function(group, data, drop=NULL, breaks=NULL, weights, grouping = NULL){
  if(!is.null(drop)){
   drop <- unique(drop)
   dropped <- match(drop, colnames(data))
   dropped <- dropped[!is.na(dropped)]

   if(length(dropped)>0) 
		data <- data[-dropped]
  }

 if(is.null(breaks)){
  vars <- colnames(data)
  nv <- length(vars)
  breaks <- vector(nv, mode="list")
  for(i in 1:nv){
   if(is.numeric(data[[i]]) | is.integer (data[[i]]))
    breaks[[i]] <- pretty(range(data[[i]],na.rm=TRUE), n=nclass.scott(na.omit(data[[i]])), 1)
   names(breaks) <- vars
  }
 }  

 if(!is.null(grouping) & !is.null(names(grouping))){
	gn <- names(grouping)
	n.gn <- length(gn)
	for(g in 1:n.gn){
		if(!is.null(data))
			data[[gn[g]]] <- group.var(data[[gn[g]]], grouping[[g]])
			if(!is.null(breaks[gn[g]]))
			breaks[gn[g]] <- NULL   			   
	}
 }
	
 n <- dim(data)[1]

 if(missing(weights)){
  weights <- rep(1, n)
 }
 rem <- which(weights<=0)
 if(length(rem)>0){
  data <- data[-rem,]
  weights <- weights[-rem]
  group <- group[-rem]
 }


 cem.imbalance(group, data, collapsed=FALSE, reduced=FALSE, breaks=breaks, weights=weights, grouping=grouping)
}

### in cem.imbalance the argument 'grouping' it is passed only as reference
### in order to be returned in L1.meas or similar. Handling of grouping
### should be done at a previous stage before calling cem.imbalance

cem.imbalance <- function(group, data, collapsed = TRUE, reduced = TRUE, breaks = NULL, 
weights, grouping = NULL)  
{    
    lv <- unique(na.omit(group))
    if (is.null(breaks)) {
        vars <- colnames(data)
        nv <- length(vars)
        breaks <- vector(nv, mode = "list")
        for (i in 1:nv) {
            if (is.numeric(data[[i]]) | is.integer(data[[i]])) 
			breaks[[i]] <- pretty(range(data[[i]], na.rm = TRUE), 
								  n = nclass.scott(na.omit(data[[i]])), 1)
            names(breaks) <- vars
        }
    }
    if (!reduced) {
        tmp <- reduce.data(data, collapse = TRUE, breaks = breaks)
        data <- tmp$data
        new.breaks <- tmp$breaks
        collapsed <- TRUE
    }
    if (!collapsed) 
	data <- collapse.data(data)
    n <- length(data)
    if (missing(weights)) {
        weights <- rep(1, n)
    }
	
    keep <- which(weights > 0)
    tmp <- data.frame(data[keep], group[keep])
    tab1 <- table(tmp)
    tmp <- unlist(apply(tab1, 1, prod))
    LCS <- length(which(tmp > 0))/length(tmp)

    nlv <- length(lv)
	for(i in 1:nlv){
		idx <- which(group == lv[i])
		weights[idx] <- weights[idx]/sum(weights[idx])
	}

    # data <- as.integer(factor(data))
	#tab <- unique(data)

    tmp <- data.frame(gr=group, w=weights,s=as.integer(factor(data)))

    #   gg <- function(x){
	#	fr <- sapply(1:nlv, function(u) sum(x$w[which(x$gr==lv[u])],na.rm=TRUE))
	#	abs(max(fr)[1] - min(fr)[1])
    #}
    #L1 <- sum(unlist(by(tmp, data, gg, simplify=TRUE)))/nlv
    
    L1 <- sum( apply(xtabs(w ~ group+s, data=tmp) , 2, function(x) abs(diff(range(x)))) )/nlv

    out <- list(L1 = L1, breaks = new.breaks, LCS = 100 * LCS, grouping=grouping)
    class(out) <- "L1.meas"
    return(out)
}


imbalance <- function(group, data, drop=NULL, breaks=NULL, weights, grouping = NULL){
 if (!is.data.frame(data))
        stop("Data must be a dataframe", call. = FALSE)

  if(!is.null(drop)){
   drop <- unique(drop)
   dropped <- match(drop, colnames(data))
   dropped <- dropped[!is.na(dropped)]

   if(length(dropped)>0) 
		data <- data[-dropped]
  }
 if(is.null(breaks)){
  vars <- colnames(data)
  nv <- length(vars)
  breaks <- vector(nv, mode="list")
  for(i in 1:nv){
   if(is.numeric(data[[i]]) | is.integer (data[[i]]))
    breaks[[i]] <- pretty(range(data[[i]],na.rm=TRUE), n=nclass.scott(na.omit(data[[i]])), 1)
   names(breaks) <- vars
  }
 }  

 if(!is.null(grouping) & !is.null(names(grouping))){
	gn <- names(grouping)
	n.gn <- length(gn)
	for(g in 1:n.gn){
		if(!is.null(data))
			data[[gn[g]]] <- group.var(data[[gn[g]]], grouping[[g]])
		if(!is.null(breaks[gn[g]]))
			breaks[gn[g]] <- NULL   			   
	}
 }
	
 n <- dim(data)[1]
 if(missing(weights)){
  weights <- rep(1, n)
 }
 rem <- which(weights<=0)
 if(length(rem)>0){
  data <- data[-rem,]
  weights <- weights[-rem]
  group <- group[-rem]
 }


 tmp <- reduce.data(data, breaks=breaks)$data
 lv <- unique(na.omit(group))
	
 globalL1 <- L1.meas(group=group, data=data, breaks=breaks, weights=weights)
	
 if(length(lv)>2){
  cat("\nMore than 2 groups, no univariate measures of imbalance will be calculated.\n")
  out <- list(L1=globalL1)
 } else {
  idx1 <- which(group==lv[1])
  idx2 <- which(group==lv[2])
 
  stdmean <- function(x, wh) weighted.mean(x, w=wh, na.rm=TRUE)

  qt.w <-function(x, probs=c(0,0.25,0.5, 0.75, 1), wh){
   q <- NULL
   which(is.na(x)) -> id
   if(length(id)>0)
    x <- x[-id]
   ord <- order(x)
   x <- x[ord]
   wh <- wh[ord]
   F <- cumsum(wh)/sum(wh)
  
   for(i in probs){
    q <- c(q, x[which(F>=i)[1]])  
   }
   idx <- which(is.na(q))
   q[idx] <- max(x)
   q
  }

  vnames <- colnames(data)
  nv <- length(vnames)
  tab <- as.data.frame(matrix(NA, nv, 8)) 
  rownames(tab) <- vnames
  colnames(tab) <- c("statistic", "type","L1", "min", "25%",  "50%",  "75%", "max")

  for (i in 1:nv){
   tab[i,3] <- L1.meas(group, tmp[i], breaks=breaks,weights=weights)$L1
   if((is.numeric(data[,i]) | is.integer(data[,i])) ){
    tab[i,1] <- stdmean(data[idx1,i],weights[idx1]) - stdmean(data[idx2,i],weights[idx2])  
    tab[i,4:8] <- qt.w(x=data[idx1,i], wh=weights[idx1]) - qt.w(x=data[idx2,i], wh=weights[idx2])  
    tab[i,2] <- "(diff)"
   } else {
    t1 <- table(data[idx1,i])
    t2 <- table(data[idx2,i])
    keep <- which(t1>0 & t2>0)
    tab[i,1] <- as.numeric(chisq.test(cbind(t1[keep],t2[keep]))$statistic)   
    tab[i,4:8] <- NA  
    tab[i,2] <- "(Chi2)"
   }
  }
  out <- list(tab=tab, L1=globalL1)
 }
 	
 class(out) <- "imbalance"
 return( out )
 
}


print.imbalance <- function(x,...){
 if(!is.null(x$L1))
  print(x$L1)
 if(!is.null(x$tab)){	
  cat("Univariate Imbalance Measures:\n\n")
  print(x$tab,...)
 }
}


print.L1.meas <- function(x,...){
 cat(sprintf("\nMultivariate Imbalance Measure: L1=%.3f", x$L1))
 cat(sprintf("\nPercentage of local common support: LCS=%.1f%%\n\n", x$LCS))
}
