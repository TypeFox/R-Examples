bw.rt<-function(x,robust=FALSE,alpha=0.5){
	x <- conversion.circular(x, units = "radians", zero = 0, rotation = "counter", modulo = "2pi")
	attr(x, "class") <- attr(x, "circularp") <- NULL
	if (!is.numeric(x)) stop("argument 'x' must be numeric")
	x.na <- is.na(x)
	if (sum(x.na)>0) warning("Missing values were removed")
	x <- x[!x.na]
	n <- length(x)
	if (n==0) stop("No observations (at least after removing missing values)")
	if (robust){
		if (!is.numeric(alpha)) stop("argument 'alpha' must be numeric")
		if (alpha<0 | alpha>1){
			warning("'alpha' must be in the interval [0,1]")
			return(NULL)
		}
		if (n<4) stop("Length of vector 'x' must be at least 4")
		x<-sort(x)
		x2<-c(x,x+2*pi)
		d<-numeric(n)
		for (i in 1:n){
			d[i]<-abs(x2[i]-x2[floor(n*alpha)+(i-1)])  
		}
		ind<-which.min(d)
		q1<-x2[ind]
		q2<-x2[floor(n*alpha)+(ind-1)] 
		q.grid<-circular(seq(as.numeric(q1),as.numeric(q2),len=500))
		mu<-(q.grid[1]+q.grid[length(q.grid)])/2
		c<-q.grid[2]-q.grid[1]
		g<- function(kappa){ alpha - (c/2 * (dvonmises(q.grid,circular(mu),kappa)[1] + dvonmises(q.grid,circular(mu),kappa)[length(q.grid)]) +
     		 c * sum(dvonmises(q.grid,circular(mu),kappa)[2:(length(q.grid)-1)])) }
		b<-10
		while(sign(g(0))== sign(g(b))){
			b<-b+10
		}
		kappa<-uniroot(g, c(0, b), tol = 0.0001)$root

	}else{
		mu <- atan2(sum(sin(x)),sum(cos(x)))
		kappa.eq<-function(kappa){besselI(kappa,1)/besselI(kappa,0)-mean(cos(x-mu))}
		b<-10
		while(sign(kappa.eq(0))== sign(kappa.eq(b))){
			b<-b+10
		}
		kappa <- uniroot(kappa.eq, c(0,b), tol = 0.0001)$root
	}
	return((3*n*kappa^2*besselI(2*kappa,2)/(4*pi^0.5*besselI(kappa,0)^2))^(2/5))
}