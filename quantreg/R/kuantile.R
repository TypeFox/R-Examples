"kuantile"<-
function(x, probs = seq(0,1, 0.25),na.rm = FALSE, 
	names = TRUE, type = 7, ...) {
    if(na.rm)
	x <- x[!is.na(x)]
    else if(any(is.na(x)))
	stop("NA's and NaN's not allowed in 'x' if 'na.rm' is FALSE")
    if(any(is.na(probs)))
	stop("NA's and NaN's in 'probs' not allowed")
    if(any(probs < 0 | probs > 1))
        stop("probs outside [0,1]")
    p <- probs
    op <- order(p)
    p <- p[op]
    n <- length(x)
    m <- length(p)
    g <- rep(.5,m)
    if(type == 1) #Hyndman-Fan Typology
	k <- j <- pmax(1,ceiling(p*n))
    else if(type == 2){
	j <- pmax(1,floor(p*n))
	k <- sort(c(pmax(1,j),pmin(j+1,n)))
	g <- ifelse(p*n > j, 1, 0.5)
	}
    else if(type == 3)
	k <- j <- pmax(1,round(p*n))
    else{
	switch(type - 3, 
		{a <- 0; b <- 1},#Type 4
		{a <- b <- 0.5}, #Type 5
		{a <- b <- 0},   #Type 6
		{a <- b <- 1},   #Type 7
		{a <- b <- 1/3}, #Type 8
		{a <- b <- 3/8}) #Type 9
	d <- a + p * (1 - a - b)
        j <- floor(p*n + d)
	g <- p*n + d - j
    	k <- sort(c(pmax(1,j),pmin(j+1,n)))
	}
    uk <- kunique(k)
    uz <- kselect(x,uk$xU)
    z <- uz[uk$ix]
    if(type %in% c(1,3))
	A <- matrix(z,m,2)
    else
	A <- t(matrix(z,2,m)) 
    G <- cbind(1-g,g)
    y <- (A * G) %*% c(1,1) # <=> diag(crossprod(A,G))
    y <- y[rank(probs)]
    if(names && m > 0){
	dig <- max(2,getOption("digits"))
	names(y) <- paste(format(100*probs, trim = TRUE, digits = dig), "%", sep="")
	}
    class(y) <- "kuantile"
    return(y)
}
"kselect" <-
function(x,k){
    n <- length(x)
    m <- length(k)
    z <- .Fortran("kuantile",
        	k = as.integer(k),
        	m = as.integer(m),
        	n = as.integer(n),
        	x = as.double(x),
		PACKAGE = "quantreg")
    return(z$x[z$k])
    }
"kunique" <-
function (x, isuniq = !duplicated(x))
{ # shamelessly plagarized from Martin Maechler's sfsmisc package
    need.sort <- is.unsorted(x)
    if (need.sort) {
        xs <- sort(x, index.return = TRUE)
        ixS <- xs$ix
        isuniq <- isuniq[ixS]
        x <- xs$x
    }
    ix <- as.integer(cumsum(isuniq))
    if (need.sort)
        ix <- ix[sort.list(ixS)]
    list(ix = ix, xU = x[isuniq])
}
