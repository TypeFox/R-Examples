

ina <- function(x){
	return(is.na(x)+0)
}

visna <- function(x, freqvar = "Freq", tp = FALSE, fr = 1, fc = 1, sort = "n", sort.method = "count",
 col = "w", mar.col = c(alpha("black",0.7),alpha("darkred",0.8),"red","green"), s = Inf, pmax = 1,  opts = list(), plot = TRUE, return.data = !plot, ...){
 	if(!any(is.na(x))){
 		stop("No NA's in the data. For indicator matrices please use visid(x, ... ) and for factor data.frames there is visdf(x,freqvar)")
 	}
 	return(
 		visid(is.na(x)+0, freqvar = freqvar, tp = tp, fr = fr, fc = fc, sort = sort, sort.method = sort.method,
 col = col, mar.col = mar.col, s = s, pmax = pmax,  opts = opts, plot = plot, return.data = return.data, ... )
 	)
 }


visdf <- function(x, freqvar = "Freq", tp = FALSE, fr = 1, fc = 1, sort = "n", sort.method = "count",
 col = "w", mar.col = c(alpha("black",0.7),alpha("darkred",0.8),"red","green"), s = Inf, pmax = 1,  opts = list(), plot = TRUE, return.data = !plot, ...){
 	if(!inherits(x,"data.frame")){
 		stop("x is not a data.frame. For indicator matrices please use visid(x, ... )")
 	}
 	return(
 		visid(idat(x, allcat = TRUE, keep = freqvar), freqvar = freqvar, tp = tp, fr = fr, fc = fc, sort = sort, sort.method = sort.method,
 col = col, mar.col = mar.col, s = s, pmax = pmax,  opts = opts, plot = plot, return.data = return.data, ... )
 	)
 }

visid <- function(x, freqvar = "Freq", tp = FALSE, fr = 1, fc = 1, sort = "n", sort.method = "count",
 col = "w", mar.col = c(alpha("black",0.7),alpha("darkred",0.8),"red","green"), s = Inf, pmax = 1,  opts = list(), plot = TRUE, return.data = !plot, ... ){
	
	if(is.null(dim(x))){
		stop("No dim attribute found: matrix or data.frame required.")
	}
	if(dim(x)[1] < 2){
		stop("x must have at least 2 rows.")
	}
	if(dim(x)[2] < 2){
		stop("x must have at least 2 columns.")
	}
	
	if("gap.prop" %in% names(opts)){
		gap.prop <- opts$gap.prop
	}else{
		gap.prop <- NULL
	}
	
	if("abbrev" %in% names(opts)){
		abbrev <- opts$abbrev
	}else{
		abbrev <- 12
	}
	
	if("bg.col" %in% names(opts)){
		bg.col <- rep(opts$bg.col,3)[1:3]
	}else{
		bg.col <- alpha(1,c(0.1,0.05,0.05))
	}
	
	if("mar" %in% names(opts)){
		mar <- opts$mar
	}else{
		mar <- 0.2
	}

	if("shape" %in% names(opts)){
		shape <- opts$shape
	}else{
		shape <- "r"
	}
	
	if("border" %in% names(opts)){
		border <- opts$border
	}else{
		border <- c(0.05,0.2)[c(tp,!tp)+1]
	}
	if(!"Freq" %in% dimnames(x)[[2]]){
		freqvar <- NULL
	}
	
# colors for margin tiles and borders

mar.col <- rep(mar.col,2)
tmp1 <- ((  col2rgb(mar.col[1]) + c(255,0,128) )/255) %% 1
tmp2 <- ((  col2rgb(mar.col[2]) + c(255,0,128) )/255) %% 1


mar.col[3] <- 2 #rgb( tmp1[1], tmp1[2] , tmp1[3]  )

mar.col[4] <- 3 #rgb( tmp2[1], tmp2[2] , tmp2[3]  )
	
	#if(is.null(freqvar)){
	#		if("Freq" %in% dimnames(x)[[2]]) freqvar <- "Freq"
	#}
	
	m <- ncol(x) - !is.null(freqvar)
	N <- nrow(x)
	
	xs <- subtable(x, 1:m, freqvar = freqvar)
	n <- nrow(xs)
	
	# get number of different category sets
		lvl <- lapply(xs[,-ncol(xs)],function(z) levels(as.factor(z)))
		mlvl <- sapply(lvl, function(z) sapply(lvl,identical,z))
		ng <- length(unique(apply(mlvl,1,paste,collapse=":")))
	
	if(is.data.frame(x)){
		x <- data.matrix(x)
	}
	
	if(is.matrix(x)){
		
		
		w <- xs[,m+1]
		dim(w) <- c(n,1)

# test code
mm <- melt(as.data.frame(xs),ncol(xs))
tt <- xtabs(mm[,1]~variable+value,data=mm)

		xs <- as.matrix(xs)[,-(m+1)]
		#xs<- as.matrix(idat( xs[,-(m+1)] , TRUE))
		
		w2 <- apply(xs,2,weighted.mean,w,na.rm=TRUE)[1:m]
		dim(w2) <- c(1,m)
		
	}

# ------------ variable groups --------------
# sep.var options
	if("var.col" %in% names(opts)){
		var.col <- opts$var.col
	}else{
		var.col <- NULL
	}
	
# get groups
	
	grps <- attr(x,"var.ids")
	if(is.null(grps)){
		#require(amap)
		# cc <- suppressWarnings(cor(xs))
		#cc[is.na(cc)] <- 0
		
		#hc <- hclust(as.dist(1+cc))
		#grps <- subtree(hc,h=1-1e-12)$data
	  grps <- 1
	}
	ngrp <- length(table(grps))
	
# ------------ variable groups --------------	

	
	
	ow <- order(w,decreasing = TRUE)
	ow2 <- order(w2,decreasing = TRUE)
	rs <- cs <- FALSE
	if( sort %in% c("b","both")){
		rs <- cs <- TRUE
	}
	if( sort %in% c("r","row","rows")){
		rs <- TRUE
	}
	if( sort %in% c("c","col","column","columns")){
		cs <- TRUE
	}
	if(sort %in% c("n","none",FALSE)){
		sort.method <- "n"
	}
	or <- 1:nrow(xs)
	oc <- 1:ncol(xs)
	
	if(sort.method %in%c( "c", "count")){
		if(rs){
			xs <- xs[ow,,drop=FALSE]
			w <- w[ow,,drop=FALSE]
			or <- ow
		}
		if(cs){
			xs <- xs[,ow2,drop=FALSE]
			w2 <- w2[,ow2,drop=FALSE]
			grps <- grps[ow2]
			oc <- ow2
		}
	}

	
	
	
	# --- filtering by nr/proportion
	
	if(fr < 1){
		fr <- max(2, min( which(cumsum(sort(w,decreasing=TRUE))/sum(w) >= fr ) ))
	}
	if(fr > 1){
		fr <- round(fr,0)
		keep <- which(rank(w,ties.method="first") > n-fr)
		xs <- xs[keep,,drop=FALSE]
		w <- w[keep,,drop=FALSE]
		or <- or[keep]
	}
	
	if(fc < 1){
		fc <- max(2, min( which(cumsum(sort(w2,decreasing=TRUE))/sum(w2) >= fc ) ))
	}
	if(fc > 1){
		fc <- round(fc,0)
		keep <- which(rank(w2,ties.method="first") > m-fc)
		xs <- xs[,keep,drop=FALSE]
		w2 <- w2[,keep,drop=FALSE]
		oc <- oc[keep]
	}
	
	if( is.null(gap.prop)){
		gap.prop <- c(0.1,0)[1 + (dim(xs) > 40)]
	}
	


if(sort.method == "ME"){
		dims <- NULL
		if(rs) dims <- 1
		if(cs) dims <- c(dims,2)
		if(!is.null(dims)){
			xs <- optME(xs,dims=dims,solver = "nn", ... )
		}
		or <- attr(xs,"orders")[[1]]
		oc <- attr(xs,"orders")[[2]]
		w <- w[or,,drop=FALSE]
		w2 <- w2[,oc,drop=FALSE]
	}

	if(sort.method == "optile"){
		perm.cat <- c(rs,cs)
		if(any(perm.cat)){
			ord <- optile(xs*w[,1],perm.cat=perm.cat,return.data=FALSE, ... )
		}
		or <- ord[[3]][[1]]
		oc <- ord[[3]][[2]]
		w <- w[or,,drop=FALSE]
		w2 <- w2[,oc,drop=FALSE]
		xs <- xs[or,oc]
	}




# scale w	
	scnd <- sort(w,decreasing=TRUE)[2]

	mark <- c(NA,mar.col[3])[ 1 + (w > s*scnd) ]
	w <- sapply(w,min,s*scnd)
	
	dim(w) <- c(length(w),1)
	
#scale w2

	mark2 <- c(NA,mar.col[4])[ 1 + (w2 > pmax) ]
	w2 <- sapply(w2,min,pmax)
	
	dim(w2) <- c(1,length(w2))


	
	# colour
	nc <- max(xs)+1
	
	# variable groups colour
if(!is.null(var.col)){
		var.colv <- getcolors(ngrp, var.col)[grps]
		tile.col <- rep(var.colv, each=nrow(xs))
		tile.col[xs == 0] <- "lightgrey"
}else{
	colv <- c("lightgrey",getcolors(nc-1,col))
	tile.col <- colv[xs+1]
	dim(tile.col) <- dim(xs)
}
	
	
	
if(plot){	
	####### PLOT #######
	grid.newpage()
	
	#rb <- lb <- bd
	border <- rep(border,2)
	mar <- rep(mar,2)
	rb <- mar[1]
	lb <- mar[2]
	
	

	if(tp){
		vm <- viewport((1-rb)/2,(1+lb)/2,1-rb,1-lb)
		pushViewport(vm)
		fluctile(t((xs+1)>0), shape = shape, add=TRUE, gap.prop = rev(gap.prop), tile.col=t(tile.col),label=c(TRUE,FALSE), border = c(border[1],0,0,border[2]),lab.opt=list(rot=0,lab.cex=0.8, abbrev = abbrev), bg.col = bg.col[1])
		popViewport()
	
	
		vr <- viewport(1-rb/2,(1+lb)/2,rb,1-lb)
		pushViewport(vr)
		fluctile(t(w2) > 0, gap.prop = rev(gap.prop), add=TRUE,dir="h",just="l", label = FALSE,tile.col=NA, border = c(0.05,0.05,0,border[2]))
		fluctile(t(w2),add=TRUE, gap.prop = rev(gap.prop), dir="h",just="l", label = FALSE, 	tile.col=mar.col[2], bg.col= bg.col[2],maxv=pmax, border = c(0.05,0.05,0,border[2]),tile.border=mark2,lab.opt=list( lwd=1+log(log(ncol(xs))) ) )
		popViewport()
	
	
		vl <- viewport((1-rb)/2,lb/2,1-rb,lb)
		pushViewport(vl)
		fluctile(t(w), gap.prop = rev(gap.prop), add=TRUE,dir="v",just="t", label = FALSE, border = c(border[1],0,0.05,0.05),tile.border=mark, tile.col = mar.col[1],lab.opt=list(lwd=1+log(log(nrow(xs)))),bg.col = bg.col[3])
		
	##vls <- viewport(0.05,0,width=0.95,height=0.9,just=c("left","bottom"))
	##rmb(tt,spine=TRUE,label=FALSE,add=TRUE,gap.prop=0.1,label.opt=list(yaxis=FALSE,s0=0),vp=vls,col=colv)
	
		popViewport()

		
	}else{
		vm <- viewport((1-rb)/2,(1+lb)/2,1-rb,1-lb)
		pushViewport(vm)
		#fluctile(xs, add=TRUE)
		fluctile((xs+1)>0, shape = shape, add=TRUE, gap.prop = gap.prop, tile.col=tile.col,label=c(FALSE,TRUE), border = c(border[1],0,0,border[2]),lab.opt=list(rot=0,lab.cex=0.8, abbrev = abbrev),bg.col = bg.col[1])
		popViewport()
	
	
		vr <- viewport(1-rb/2,(1+lb)/2,rb,1-lb)
		pushViewport(vr)
		fluctile(w, gap.prop = gap.prop, add=TRUE,dir="h",just="l", label = FALSE, border = c(0.05,0.05,0,border[2]),tile.border=mark, tile.col = mar.col[1],lab.opt=list(lwd=1+log(log(nrow(xs)))),bg.col = bg.col[3])
		popViewport()
	
	
		vl <- viewport((1-rb)/2,lb/2,1-rb,lb)
		pushViewport(vl)
		fluctile(w2 > 0, gap.prop = gap.prop, add=TRUE,dir="v",just="t", label = FALSE,tile.col=NA, border = c(border[1],0,0.05,0.05))
		fluctile(w2,add=TRUE, gap.prop = gap.prop, dir="v",just="t", label = FALSE, tile.col=mar.col[2],bg.col = bg.col[2],maxv=pmax, border = c(border[1],0,0.05,0.05),tile.border=mark2,lab.opt=list(lwd=1+log(log(ncol(xs)))))
	##vls <- viewport(0.05,0,width=0.95,height=0.9,just=c("left","bottom"))
	##rmb(tt,spine=TRUE,label=FALSE,add=TRUE,gap.prop=0.1,label.opt=list(yaxis=FALSE,s0=0),vp=vls,col=colv)
	
		popViewport()

		
		
	}
} # END IF plot
if(return.data){
	attr(xs,"mar") <- list(rm = w, cm = w2)
	attr(xs,"orders") <- list(or,oc)
	return(xs)
}	
	
	return(invisible(TRUE))
}



ibmat <- function(x, numeric = 0.2,  df = 5){
	#require(splines)
	if( !( "bs" %in% ls() ) ) bs <- function(x) x^0
	
		if(numeric[1] < 1){
 			numeric <- sapply(x,function(z){
 				if(is.factor(z)){
 					return(FALSE)
 				} 
 				ret <- any( abs(z)%%1 > 0) || length(table(z)) > numeric*nrow(x)
 				return(ret)
 			})
 		}
 			
 		D1 <- idat(x[,!numeric],allcat=TRUE)
 			
		D2 <- do.call(cbind,lapply(x[,numeric], bs, df = df))
 			
		if(length(D1) == 0) D1 <- NULL
		if(length(D2) == 0){
			D2 <- NULL
		}else{
			names(D2) <- paste(1:df,rep(names(x)[numeric],each=df),sep=":")
		} 
		
		
		D <- data.frame(D1,D2)
	
	return(D)
}


