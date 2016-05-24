#Functions used internally by more than one metRology function

.get.pars<-function(distrib.name, x, u, df=NULL) {

	supported.distribs<-c("norm", "unif", "tri", "t", "t.scaled")
	
	if(!(distrib.name %in% supported.distribs) ) {
		stop(paste("Unable to get parameters for:", distrib.name), call. = TRUE)
	}
	
	if( (distrib.name %in% c("t", "t.scaled")) && !is.numeric(df) ) {
		stop(paste(distrib.name, "distribution requires numeric df", sep=" "), call. = TRUE)
	}
	
	switch( distrib.name,
		norm=list(mean=x, sd=u),
		unif=list(min=x-sqrt(3)*u, max=x+sqrt(3)*u),
		tri=list(min=x-sqrt(6)*u, max=x+sqrt(6)*u, mode=x),
		t=list(df=df, mean=x, sd=u),
		t.scaled=list(df=df, mean=x, sd=u)
	)
	
}



.get.df<-function(k, p) {
	cl <- 1-(1-p)/2
	if( k <= qnorm(cl) ) {
		return(Inf)
	} else {
		d <- 1
		i <- 1
		while(qt(cl, d+i) > k ) {
			i <- 2*i
		}
		d<-round(uniroot(function(d, k) qt(cl, d) - k , c(d, d+i), k=k)$root, 0)
		return(d)
	}
}

