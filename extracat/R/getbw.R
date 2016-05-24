getbw <- function(x, k = NULL, min_n = NULL, warn = FALSE){

N <- length(x)
	if(is.null(k)) k <- 1 + 2*ceiling(log(N)/log(2))
	
	if(is.null(min_n)) min_n <- max(log(N/10)/log(10),1)

if( (min_n > 0) & (k*min_n >  length(x)/min_n)){
	rat <- k*min_n / length(x)
	print(rat)
	a <- round(rat*min_n)
	if( a == 0 ) a <- 1
	min_n <- a
	if(warn) cat("The number of observations is smaller than k*min_n. Reducing min_n to",min_n) 	
}
if(k > (k2<-length(unique(x))) & min_n > 0){
	k <- k2
	if(warn) cat("There are too few different values in x. Setting k to ",k2) 	
}	

rg <- range(x)

x <- sort(x)
dx <- diff(x)
ut <- min(dx[which(dx>0)])


span <- diff(rg)
bw <- span/ut/k
old.bw <- bw
anchor <- min(x)

optimal <- FALSE
while(!optimal){
	
	cx <- cut(x, breaks <- seq(anchor, max(x)+bw*ut,bw*ut),right = FALSE)
	tx <- table(cx, useNA = "ifany")
	curr.k <- length(which(tx >= min_n))
	#if( curr.k < k ){
			mult <- curr.k/k
			old.bw <- bw
			bw <- bw*mult
			
	#}else{
		if(curr.k == k){
			optimal <- TRUE
		}#else{
		#	bw <- round((old.bw+bw)/2,0)
		#}
			
	#}
		if(mult < 1 & bw < 1){
			optimal <- TRUE
		}
		
	}
bwf <- max(1,floor(bw))

bwc <- ceiling(bw)
	bw <- bwc
	cx1 <- cut(x, breaks <- seq(anchor, max(x)+bw*ut,bw*ut),right = FALSE)
	tx1 <- table(cx1, useNA = "ifany")
	curr.k1 <- length(which(tx1 >= min_n))
	tx <- tx1
	if(curr.k1 != k){
		bw <- bwf
		cx <- cut(x, breaks <- seq(anchor, max(x)+bw*ut,bw*ut),right = FALSE)
		tx <- table(cx, useNA = "ifany")
		curr.k <- length(which(tx >= min_n))
		if(abs(curr.k-k) > abs(curr.k1-k)){
			curr.k <- curr.k1
			bw <- bwc	
		} 
	}
	if(warn){
		if(curr.k != k){
			cat("Best solution found has k = ",curr.k, " bins.")	
		}
	}
ret <- seq(anchor,max(x)+bw*ut,bw*ut)
attr(ret,"bw") <- bw*ut
attr(ret,"k") <- curr.k
attr(ret,"outlier") <- tx > 0 & tx < min_n
return(ret)


	

}


cutbw <- function(x, k = NULL, min_n = NULL,  warn = FALSE){
	return( cut(x,	breaks=getbw(x,k,min_n,warn)) )
}


ahist <- function(x, k = NULL, m = NULL, fun = "qplot", col = "grey", ival = NULL){
	
	if(!is.null(ival)){
		if(ival == TRUE) ival <- 0.95
		stopifnot( 0 < ival & ival < 1)
		ii <- innerval(x,p=ival)
		x <- x[which(x >= ii[1] & x <= ii[2])]
	}
	
	br <- getbw(x, k = k, min_n = m)
	if(fun=="qplot"){
			p1 <- ggplot(data = data.frame(x), aes(x=x))
			p1 <- p1 + geom_bar(colour=col,binwidth = attr(br,"bw"))
			
		return(p1)
	}
	if(fun=="hist"){
			hist(x, breaks = br, col = col)
		
	}
	return(invisible(br))
}
