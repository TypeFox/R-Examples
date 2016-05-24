#Implements DerSimonian-Laird estimate
#This implementation parallels mpaule in metRology
#in allowing a vector of raw data, with grouping factor

dsl<-function(x, ..., na.rm=FALSE) {
	UseMethod("dsl")
}

dsl.default<-function(x, s, n=NULL, groups=NULL, ..., na.rm=FALSE) {
	#x is a vector of (mean) observations
	#s is a vector of standard errors or
	#(if n is given) standard deviations
	#n is an optional mumber of observations per group
	#groups is a grouping factor, used to group x
	
	#This uses the implementation given by
	#Jackson et al, J Stat Plan Inf 140, 961-970 (2010)
	
	#Prepare to remove NA's 
	if(!is.null(n)) {
		n.in <- n
		n <- rep(n, length(x)) #Recycled if present
	}
	x.in <- x 
	s.in <- if(missing(s)) NA else s
	g.in <- if(is.null(groups)) NA else groups
	
	if(na.rm) {
		na.xs <- if(missing(s)) 
				is.na(x) 
			else 
				(is.na(x) | is.na(s))
		x <- x[!na.xs]
		if(!missing(s)) s <- s[!na.xs]
		if(!is.null(groups)) groups <- groups[!na.xs]
		if(!is.null(n)) n <- n[!na.xs]
	}
	
	if(!is.null(groups)) {
		if(!missing(s) || !is.null(n)) 
			warning("Using groups to calculate stderr: ignoring s and n")

		count<-function(x) length( x )
		groups <- factor(groups)
		x <- as.vector(tapply(xi<-x, groups, mean, na.rm=na.rm))
		s <- as.vector(tapply(xi, groups, sd, na.rm=na.rm))
		n <- as.vector(tapply(xi, groups, count))
		std.err <- s/sqrt(n)
	} else {
		if(missing(s)) stop("Either s or groups must be specified")
		std.err <- if(!is.null(n)) s/sqrt(n) else s
	}
	
	
	
	N <- length(x)
	
	w <- 1/std.err^2
	
	x.bar <- sum(w*x)/sum(w)
	
	Q <- sum(w * (x - x.bar)^2)

	S1 <- sum(w)
	S2 <- sum(w^2)
	
	tau.DL.2 <- max(0, (Q-(N-1))/(S1 - S2/S1))
	
	s.adj <- sqrt(std.err^2 + tau.DL.2)
	
	mu.DL <- sum(x/s.adj^2) /
 			sum(1/s.adj^2)
 	
	s.DL <- sqrt(1/sum(1/s.adj^2))
	

	rv <- .construct.loc.est( x=mu.DL, u=s.DL, df=length(x)-1, 
		xi=x, ui=s, u.eff=s.adj, w=rep(1, length(x)), method="DerSimonian-Laird", 
		method.details=list(mu=mu.DL, s=s.DL, tau=sqrt(tau.DL.2)) )
	
	return(rv)
	
}

