bootstrapCapwire <-
function(x, bootstraps=1000, CI=c(0.025, 0.975)){
	
	if (x$model == "Equal.capture"){
		
		n <- x$ml.pop.size
		s <- x$sample.size
		m <- x$max.pop
		c <- x$cap.ind
		t <- x$sampled.ind
		
		boot.ecm <- function(n, s, m){
			x <- simEcm(n, s)
			
			if (all(unique(x[,2]) == 1)) {
				y <- m
			} else {	
			y <- suppressWarnings(fitEcm(x, m)$ml.pop.size)
			}	
			return(y)
		}
		
		boot <- sapply(1:bootstraps, function(x) boot.ecm(n, s, m))
		conf.int <- quantile(boot, CI)
		
		for (i in 1:length(conf.int)){
			if (conf.int[i] < t ){
				conf.int[i] <- t
			}
		}
		
		ml.pop.size <- n
		
		return(list(ml.pop.size=ml.pop.size, conf.int=conf.int))
	}
	
	if (x$model == "Two.innate.rates"){
		
		na <- x$ml.na
		nb <- x$ml.nb
		s <- x$sample.size
		a <- x$alpha
		m <- x$max.pop
		c <- x$cap.ind
		t <- x$sampled.ind
		
		boot.tirm <- function(na, nb, a, s, m){
			x <- simTirm(na, nb, a, s)
			
			if (all(unique(x[,2]) == 1)) {
				y <- m
			} else {	
			y <- suppressWarnings(fitTirm(x, m)$ml.pop.size)
			}
			return(y)
		}
		
		boot <- sapply(1:bootstraps, function(x) boot.tirm(na, nb, a, s, m))
		conf.int <- quantile(boot, CI)
		
		for (i in 1:length(conf.int)){
			if (conf.int[i] < t ){
				conf.int[i] <- t
			}
		}
		
		
		ml.pop.size <- x$ml.pop.size
		
		return(list(ml.pop.size=ml.pop.size, conf.int=conf.int))
	}
	
	if (x$model == "Two.innate.rates.partitioned"){
		
		na <- x$ml.na
		nb <- x$ml.nb
		a <- x$alpha
		s <- x$sample.size
		m <- x$max.pop
		c <- x$cap.ind
		t <- x$sampled.ind
		excluded <- x$excluded

		ex <- sum(excluded[,2])

		boot.tirm <- function(na, nb, a, s, m){
			x <- simTirm(na, nb, a, s)
			
			if (all(unique(x[,2]) == 1)) {
				y <- m
			} else {	
			y <- suppressWarnings(fitTirm(x, m)$ml.pop.size)
			}
			return(y)
		}
		
		
		
		boot <- sapply(1:bootstraps, function(x) boot.tirm(na, nb, a, s, m))
		conf.int <- quantile(boot, CI) + ex
		
		for (i in 1:length(conf.int)){
			if (conf.int[i] < t ){
				conf.int[i] <- t
			}
		}
		
		
		ml.pop.size <- x$ml.pop.size
		
		return(list(ml.pop.size=ml.pop.size, conf.int=conf.int))

	}
}
