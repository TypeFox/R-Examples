boot.ratio.test <- function(boot.cube,critical.value=2){	
	boot.cube.mean <- apply(boot.cube,c(1,2),mean)
	boot.cube.mean_repeat <- array(boot.cube.mean,dim=c(dim(boot.cube)))
	boot.cube.dev <- (boot.cube - boot.cube.mean_repeat)^2
	#s.boot<-(apply(boot.cube.dev,c(1,2),mean))^(1/2)
	s.boot<-sqrt(apply(boot.cube.dev,c(1,2),mean))
	boot.ratios <- boot.cube.mean / s.boot
	boot.ratios <-replace(boot.ratios,is.infinite(boot.ratios),0)
	significant.boot.ratios <- (abs(boot.ratios) > critical.value)
	rownames(boot.ratios) <- rownames(boot.cube)
	rownames(significant.boot.ratios) <- rownames(boot.cube)	
	return(list(sig.boot.ratios=significant.boot.ratios,boot.ratios=boot.ratios,critical.value=critical.value))
}