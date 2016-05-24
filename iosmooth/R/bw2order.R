bw2order <-
function(acf.taper, n) { #Plug in bandwidth for Parzen kernel estimator
	#acf.taper is the tapered acf given by the infinite order estimate
	#n is the sample size.
	int.f.sq <- (2*sum(acf.taper^2)-acf.taper[1]^2)/(2*pi) #flat top estimat of \int f(\omega)^2 d\omega
	int.lambda.sq <- 151/280 #Parzen kernel square integral  By Mathematica.
	
	int.f2d.sq <- sum(((0:(length(acf.taper)-1))^4)*acf.taper^2)/pi  #flat top estimate of \int f''(\omega)^2
	Lambda.2nd.moment <- 12 #Parzen kernel second moment.  By Mathematica.
	
	CV <- int.f.sq * int.lambda.sq
	CB <- int.f2d.sq*Lambda.2nd.moment^2/4  #Coef on h^4 in IMSE expression
	
	h <- (CV/(4*CB))^(1/5) * n^(-1/5)
	M <- ceiling(1/h)
	M
}
