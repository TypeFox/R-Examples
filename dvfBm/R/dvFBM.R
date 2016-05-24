


filt<-function(nm="i2"){
	## gives the coefficients of an increments type filter, a Daublet, Symmlet or Coiflet filter
	fact<-function(n) ifelse(n==0,1,prod(1:n))
	cnk<-function(n,k) fact(n)/fact(k)/fact(n-k)
	
	if (!is.na(o<-as.numeric(strsplit(nm,"i")[[1]][2]))) {
		l<-o+1
		a<-rep(0,l)
		for (k in 0:o){
			a[k+1]<-cnk(o,k)*(-1)^(k)
		}
	}
	else {
		require(wmtsa)
		a<-wavDaubechies(nm)$wav
	}
	a
}

dilatation <- function(a=c(1,-2,1), m=2){
	## provides the dilated version of a filter a
	if(m>1){
		la <- length(a)
        	am <- rep(0, m * la - 1)
        	am[seq(1, m * la - 1, by = m)] <- a
        	am 
	} else a
}



dvFBM<-function(fbm,nma="i2",M1=1,M2=5,method=c("ST","Q","TM","B1-ST","B1-Q","B1-TM","B0-ST","B0-Q","B0-TM"),par=list(),llplot=FALSE){

	if (missing(fbm)) stop("Missing data")
	
	listDaub<-c("d2", "d4", "d6", "d8","d10", "d12", "d14", "d16", "d18", "d20","s2","s4", "s6", "s8", "s10","s12", "s14", "s16", "s18", "s20", "l2","l4", "l6", "l14", "l18", "l20","c6", "c12", "c18", "c24", "c30")
	l<-strsplit(nma,"")[[1]][1]
	if (l=="i") {
		if (length(strsplit(nma,"i")[[1]])>2) stop("Bad entry name for the filter ")
		if (is.na(as.numeric(strsplit(nma,"i")[[1]][2]))) stop("Bad entry name for the filter")
		a<-filt(nma)
	}else{
		if (nma %in% listDaub) {
			a<-filt(nma)
		}else stop("Bad entry name for the filter")
	}

	if (method %in% c("Q","B0-Q","B1-Q")){
		if (is.null(par$vecp) | is.null(par$vecc)) stop('"par=list(vecp=...,vecc=...) is needed for methods "Q","B0-Q","B1-Q"')
	}
	if (method %in% c("TM","B0-TM","B1-TM")){
		if (is.null(par$beta1) | is.null(par$beta2)) stop('par=list(beta1=...,beta2=...) is needed for methods "TM","B0-TM","B1-TM"')
	}

	if (!(method %in% c("ST","Q","TM","B1-ST","B1-Q","B1-TM","B0-ST","B0-Q","B0-TM"))) stop('Method should be one of "ST","Q","TM","B1-ST","B1-Q","B1-TM","B0-ST","B0-Q","B0-TM"')
	
	
	l<-length(a)-1
	n<-length(fbm)
	Unam<-NULL
	Unam<-switch(method,
		"ST"={
			for (m in M1:M2){
 				am<-dilatation(a,m)
				Vam <- filter(fbm, am, sides = 1)
        			Vam <- Vam[!is.na(Vam)]
				Unam<-c(Unam,mean(Vam^2))
			}
			Unam
			
		},
		"Q"={
			vecp<-par$vecp
			vecc<-par$vecc
			for (m in M1:M2) {
      				am<-dilatation(a,m)
				Vam <- filter(fbm,dilatation(a,m),sides=1)
      				Vam <- Vam[!is.na(Vam)]
      				Unam <- c(Unam, sum(vecc * quantile(Vam^2,vecp)))
			}
			Unam
		},
		"TM"={
			beta1<-par$beta1
			beta2<-par$beta2
			for (m in M1:M2) {
			      am<-dilatation(a,m)
			      Vam <- filter(fbm, am, sides = 1)
        		      Vam <- Vam[ - (1:(m*l))]
			      tmp<-sort(Vam^2)
      			      nn<-length(Vam)
			      Unam<-c(Unam,mean(tmp[(trunc(nn*beta1)+1):(nn-trunc(nn*beta2))]))
    			}
			Unam
		},
		"B1-ST"={
			for (m in M1:M2) {
			    Va2m <- filter(fbm,dilatation(a,2*m),sides=1)
			    Va2m <- Va2m[!is.na(Va2m)]
			    Vam <- filter(fbm,dilatation(a,m),sides=1)
			    Vam <- Vam[!is.na(Vam)]
			    Unam <- c(Unam,abs(mean(Va2m^2)-mean(Vam^2)))
			}
			Unam
	    	},
		"B1-Q"={
			vecp<-par$vecp
			vecc<-par$vecc
			for (m in M1:M2) {
      				Va2m <- filter(fbm,dilatation(a,2*m),sides=1)
	        		Va2m <- Va2m[!is.na(Va2m)]
			        Vam <- filter(fbm,dilatation(a,m),sides=1)
			        Vam <- Vam[!is.na(Vam)]
      				Unam <- c(Unam, abs(sum(vecc * quantile(Va2m^2,vecp)) - sum(vecc * quantile(Vam^2,vecp))))
			}
			Unam
    		},
		"B1-TM"={
			beta1<-par$beta1
			beta2<-par$beta2
			for (m in M1:M2) {
			      Va2m <- filter(fbm,dilatation(a,2*m),sides=1)
	                      Va2m <- Va2m[!is.na(Va2m)]
			      Vam <- filter(fbm,dilatation(a,m),sides=1)
			      Vam <- Vam[!is.na(Vam)]
			      tmp<-sort(Vam^2);nn<-length(Vam)
			      tmp2<-sort(Va2m^2);nn2<-length(Va2m)
			      Unam<-c(Unam, abs(mean(tmp2[(trunc(nn2*beta1)+1):(nn2-trunc(nn2*beta2))]) -mean(tmp[(trunc(nn*beta1)+1):(nn-trunc(nn*beta2))])))
    			}
			Unam
		},
		"B0-ST"={
			for (m in M1:M2) {
			    Va2m <- filter(fbm,dilatation(a,2*m),sides=1)
			    Va2m <- Va2m[!is.na(Va2m)]
			    Vam <- filter(fbm,dilatation(a,m),sides=1)
			    Vam <- Vam[!is.na(Vam)]
			    Unam <- c(Unam,abs(mean(Va2m^2)/(2*m)-mean(Vam^2)/m))
			}
			Unam
	    	},
		"B0-Q"={
			vecp<-par$vecp
			vecc<-par$vecc
			for (m in M1:M2) {
      				Va2m <- filter(fbm,dilatation(a,2*m),sides=1)
	        		Va2m <- Va2m[!is.na(Va2m)]
			        Vam <- filter(fbm,dilatation(a,m),sides=1)
			        Vam <- Vam[!is.na(Vam)]
      				Unam <- c(Unam, abs(sum(vecc * quantile(Va2m^2,vecp))/(2*m) - sum(vecc * quantile(Vam^2,vecp))/m))
			}
			Unam
    		},
		"B0-TM"={
			beta1<-par$beta1
			beta2<-par$beta2
			for (m in M1:M2) {
			      Va2m <- filter(fbm,dilatation(a,2*m),sides=1)
	                      Va2m <- Va2m[!is.na(Va2m)]
			      Vam <- filter(fbm,dilatation(a,m),sides=1)
			      Vam <- Vam[!is.na(Vam)]
			      tmp<-sort(Vam^2);nn<-length(Vam)
			      tmp2<-sort(Va2m^2);nn2<-length(Va2m)
			      Unam<-c(Unam, abs(mean(tmp2[(trunc(nn2*beta1)+1):(nn2-trunc(nn2*beta2))])/(2*m) -mean(tmp[(trunc(nn*beta1)+1):(nn-trunc(nn*beta2))])/m ))
    			}
			Unam
		}
)
	reg<-lm(log(Unam)~log(M1:M2))
	opt<-rev(reg$coeff)
	
	if (method %in% c("B0-ST","B0-Q","B0-TM")) Hest<-(opt[1]+1)/2
	else Hest<-opt[1]/2

	if (llplot) { plot(log(Unam)~log(M1:M2));abline(reg)}
	Hest
}
	
