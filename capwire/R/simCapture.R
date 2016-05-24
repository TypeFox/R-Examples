simCapture <-
function(n, s, dist.func, return.cap.probs=FALSE){
	
	sp <- dist.func(n)
	
	dd <- sample(c(1:n), s, replace=TRUE, prob=sp) ## sample populations
	
	get.cap <- function(ind, cap.dat, cap.prob){
		
		x <- length(cap.dat[cap.dat == ind])
		out <- cbind(ind, x, cap.prob)
		return(out)		
	}
	
	dat <- lapply(c(1:n), function(x) get.cap(x, dd, sp[x]))
	dat <- do.call(rbind, dat)
	colnames(dat) <- c("ID", "No.Obs", "Cap.Prob")
	
	if (return.cap.probs == TRUE){
		
		ind.data <- dat
		
		class.data <- indToClass(dat)
		
		return(list(ind.data=ind.data, class.data=class.data))
	
	}
	
	dat <- dat[-which(dat[,2] == 0),c(1,2)]
		
	data <- indToClass(dat)
	
	return(data)
	
}
