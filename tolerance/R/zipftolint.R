zipftol.int <- function (x, m = NULL, N = NULL, alpha = 0.05, P = 0.99, side = 1, 
	s = 1, b = 1, dist = c("Zipf", "Zipf-Man", "Zeta"), ties = FALSE, ...) {
    dist = match.arg(dist)
    if (side != 1 && side != 2) {
        stop(paste("Must specify a one-sided or two-sided procedure!", 
            "\n"))
    }
    if (side == 2) {
        alpha <- alpha/2
        P <- (P + 1)/2
    }
	fit <- zm.ll(x = x, N = N, s = s, b = b, dist = dist, ...)
	if(class(x)!="table"){ 
	x=table(x)
	}
	names(x)=1:length(x)
	x.labs=as.numeric(names(x))
	N.temp=max(x.labs)
	if(is.null(N)){
	N=N.temp
	}
	n <- sum(x)
	if(is.null(m)) m<-n
	if(dist=="Zipf"){	
		s.hat <- as.numeric(stats4::coef(fit))
		s.se <- sqrt(as.numeric(stats4::vcov(fit)))
		CI <- s.hat + c(-1,1)*qnorm(1-alpha)*s.se*sqrt(n/m)
		lower.s <- max(CI[1],0)
		upper.s <- CI[2]
		lower <- max(qzipfman(1 - P, s = upper.s, N = N), 1)
    	upper <- min(qzipfman(P, s = lower.s, N = N), N)
	} 
	if(dist=="Zipf-Man"){
		s.hat <- as.numeric(stats4::coef(fit)[1])
		s.se <- sqrt(as.numeric(stats4::vcov(fit)[1,1]))
		b.hat <- as.numeric(stats4::coef(fit)[2])
		b.se <- sqrt(as.numeric(stats4::vcov(fit)[2,2]))
    		if (b.hat==0) warning("MLE for b is 0! Consider fitting a Zipf distribution.")
		s.CI <- s.hat + c(-1,1)*qnorm(1-alpha)*s.se*sqrt(n/m) #These intervals might need to be divided by 2 for joint confidence.
		b.CI <- b.hat + c(-1,1)*qnorm(1-alpha)*b.se*sqrt(n/m)
		lower.s <- max(s.CI[1],1e-14)
		upper.s <- s.CI[2]
		lower.b <- max(b.CI[1],0)
		upper.b <- b.CI[2]	
    		lower <- max(qzipfman(1 - P, s = upper.s, b = lower.b, N = N), 1)
    		upper <- min(qzipfman(P, s = lower.s, b = upper.b, N = N), N)
	}
	if(dist=="Zeta"){
		N <- Inf
		s.hat <- as.numeric(stats4::coef(fit))
		s.se <- sqrt(as.numeric(stats4::vcov(fit)))
		CI <- s.hat + c(-1,1)*qnorm(1-alpha)*s.se*sqrt(n/m)
		lower.s <- max(CI[1],0)
		upper.s <- CI[2]
	    	lower <- max(qzipfman(1 - P, s = upper.s, N = Inf), 1)
    		upper <- qzipfman(P, s = lower.s, N = Inf)
	}	
    	if (class(lower) == "try-error") {
	        lower <- 1
	}    
	if (class(upper) == "try-error") {
	        upper <- N
	}
	if (side == 2) {
	        alpha <- 2 * alpha
	        P <- (2 * P) - 1
	}
	if(ties & length(x) > upper){
	        upper <- max(as.numeric(which(x==x[upper])))
	}
	if(ties & lower > 1){
	        lower <- min(as.numeric(which(x==x[lower])))
	}
	if(dist=="Zipf-Man"){
		temp <- data.frame(cbind(alpha, P, s.hat, b.hat, lower, upper))
		if (side == 2) {
		        colnames(temp) <- c("alpha", "P", "s.hat", "b.hat", "2-sided.lower","2-sided.upper")
	} else {
		colnames(temp) <- c("alpha", "P", "s.hat", "b.hat", "1-sided.lower", "1-sided.upper")
		}    
	} else{
		temp <- data.frame(cbind(alpha, P, s.hat, lower, upper))
		if (side == 2) {
			colnames(temp) <- c("alpha", "P", "s.hat", "2-sided.lower","2-sided.upper")
		} else {
			colnames(temp) <- c("alpha", "P", "s.hat", "1-sided.lower","1-sided.upper")
		}
	}
temp
}