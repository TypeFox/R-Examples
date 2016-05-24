.packageName <- "kza"

kz <- function(x,m,k = 3) {
    if (length(dim(x)) > 3) stop("Too many dimensions.")
    if (is.ts(x)) {
    	TS=TRUE
    	start=start(x)
    	f=frequency(x)
    } else {TS=FALSE}
    storage.mode(x) <- "double"
    x <- .Call("kz", x, as.vector(as.integer(floor(m/2))), as.integer(k))
    if (TS) {
    	x<-ts(x, start=start, frequency=f)
    }
    return (x)
}

kza <- function(x, m, y = NULL, k = 3, min_size = round(0.05*m), tol = 1.0e-5, impute_tails = FALSE) {
    if (length(dim(x)) > 3) stop("Too many dimensions.")
    if (is.null(y)) y<-kz(x,m=m,k=k)

    if (is.ts(x)) {
    	TS=TRUE
    	start=start(x)
    	f=frequency(x)
    } else {TS=FALSE}
    
    storage.mode(x) <- storage.mode(y) <- "double"
    kza.x <- .Call("kza", x, y, as.vector(as.integer(floor(m))), as.integer(k), as.integer(min_size), as.double(tol))
    if (TS) {
    	kza.x<-ts(kza.x, start=start, frequency=f)
    }
    
    if (impute_tails==FALSE) {
    	kza.x[1:m]=NA
    	kza.x[(length(kza.x)-m):length(kza.x)]=NA
	}    	
    
 	structure(list(
 			time.series = x,
            kz = y,
            kza = kza.x,
            window=m, k=k, min_size=min_size, tol=tol,
            call=match.call()
            ),
        class = "kza")    
}

plot.kza <- function(x, ...)
{
	if (is.ts(x$kz) && is.ts(x$kza)) {
    	plot(cbind(
               kz = x$kz,
               kza = x$kza
               ),
         main = paste("KZA Decomposition of time series"), 
         ...)
	} else {
		par(mfrow=c(2,1))
		plot(x$kz, ylab="kz", type='l')
		plot(x$kza, ylab="kza", type='l')
		par(mfrow=c(1,1))
	}         
}

kzsv <- function(object) 
{
	if (class(object) == "kza") {y=object}
	else { stop("Need to use result from kza!") }
	
    s <- .Call("kzsv", y$kza, y$kz, y$window, y$k, y$m, y$tol)

    if (is.ts(y$kz)) s<-ts(s, frequency=frequency(y$kz), start=start(y$kz))

 	structure(list(
 			kza = y,
            kzsv = s,
            call=match.call()
            ),
        class = "kzsv")    
}

plot.kzp <- function(x, ...)
{
	if (is.null(x$smooth_periodogram)) dz<-x$periodogram else dz<-x$smooth_periodogram
	omega<-seq(0:(length(x$periodogram)-1))/(x$k*(x$window-1))
	plot(omega, dz, type="l", xlab="Frequency", ylab="")
}

plot.kzsv <- function(x, ...)
{
	x<-x$kza
  plot(cbind(	kz = x$kz,	kza = x$kza, sigma= sqrt(x$kzsv/mean(x$kzsv))/1.96, 
    concavity=diff(diff(sqrt(x$kzsv/mean(x$kzsv))/1.96))), type='l',	    	
   		main = paste("KZSV Sample Variance"))
}

.cluster <- function(x, span=15)
{
	p=NULL
	m=NULL
	i=1
	z=x[i]
	for (y in x) {
		if (abs(y-z)<span) p=c(p,y) else  {m=c(m,round(mean(p))); p=c(y); } 
		z=x[i]
		i=i+1
	}
	m=c(m,round(mean(p)))
}

.peaks <- function(x, sigma=3, span=25)
{
	a<-x-mean(x)
	p2=NULL
	i=0

	sigma=sigma*sqrt(var(x))
	p<-a[a>sigma]
	p2<-rep(0,length(p))
	i=1
	for (j in p) {
		p2[i] <-	which(a==j)
		i=i+1
	}
	return (.cluster(p2, span))
}


summary.kzsv <- function(object, digits = getOption("digits"), ...)
{
    cat(" Call:\n ")
    dput(object$call, control=NULL)

	s<-sqrt(object$kzsv/mean(object$kzsv))/1.96
	d<-diff(diff(sqrt(object$kzsv/mean(object$kzsv))/1.96))
	
	p<-.peaks(s)
	
	if (is.ts(object$kzsv)) {
	    cat("\n Dates of interest:\n")
	    cat(" dates \t\t sigma\n")
	    for (m in p) {
	    	cat (" "); cat(as.integer(time(s)[m])); cat("\t");
	    	cat(" "); cat(cycle(s)[m]); cat("\t");
	    	cat(" "); cat(round(s[m],1)); cat("\n");
	    }
	} else {
	    cat("\n Periods of interest:\n")
	    cat(" period\n")
	    cat(" "); cat(p); cat("\t\t"); cat("\n");
	}
	
    invisible(object)
}

kzs <- function(y,m=NULL,k=3,t=NULL) 
{
	if (is.null(m)) {
		m=100*sqrt(mean(diff(y)^2))
		if (m>length(y)) m=2
	}		
	
	y<-kzft(y,m=m,k=k,f=0)
	a<-y

	if (is.ts(y)) ans<-ts(Re(a),start=start(y),frequency=frequency(y))
	else ans<-Re(a)
	return (ans)
}

