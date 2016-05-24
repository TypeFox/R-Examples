eseq <- function(n){
	if (n %% 2 != 0){
	   e <- 1
	   for (i in seq(3,to=n,by=2)){
	   alt <-rep((i-1):i,length.out=i-2)
	   new <- as.vector(rbind(1:(i-2),alt))
	   new <- c(new,i,1)
	   e <- c(e,new[-1])
	   }}
	else {
	e <- 1:2
	   if (n >= 4)
	   for (i in seq(4,to=n,by=2)){
	   alt <-rep((i-1):i,length.out=i-2)
	   new <- as.vector(rbind(alt,1:(i-2)))
	   new <- c(new,i-1,i)
	   e <- c(e,new)
		}}
	return(e)
	}
	

	   	
eseqa <- function(n){
	if (n %% 2 != 0) {
	  	m <- (n-1)/2
	    e <- array(0,c(n,m))
	    e[1,1] <- 0
	    for (k in 2:n)  
		    e[k,1] <- e[k-1,1] + m*(m+1)/2  
        if (m >= 2) for (j in 2:m) e[,j] <- e[,j-1] + j-1
        e <- e%% n +1	
        e <- c(as.vector(t(e)),1)
	  	}
	else {
		e <- eseqa(n+1)
		e <- e[e!=(n+1)]
		e <- e[-length(e)]
		}
	return(e)
	}



kntour_drop <- function(e){
	n <- max(e)
	if (n %% 2 != 0){
	# drop n from an euler seq for odd n
    m <- length(e)
	
	if (n==e[1])  x <- e[e!=n] 
	else if (n==e[m-1]) {
		x <- e[-m]
		x <- x[x!=n]}
	else if (n==e[2]) {
		x <- e[-1]
		x <- x[x!=n] }
	else {
	    x<- e[-m]
	    m <- length(x)
	    i <- match(n,x)	    
	    x <- x[c(i:m,1:(i-1))]
	    x <- x[x!=n]}
	return(x)}
	else {
	
	stop("Argument must be an euler tour on 1..n for odd n")}
	}
	
kntour_add <- function(e){
	n <- max(e)
	if (n %% 2 != 0){
	# add n+1 to an euler seq for odd n
     new <- 1:n
     new <- new[new!=e[1]]
     enew <- n+1
     for (j in seq(2,n-1,2)) enew <- c(enew,new[(j-1):j],n+1)
     return(c(e,enew))}

	else {
	
	stop("Argument must be an euler tour on 1..n for odd n")}
}
