findgacv.scad <- function(y,model){
	#
	# By Axel Benner (15.03.2007)
	
	#make shure that y is a vector of numeric numbers! 
	# and have values -1 and 1
	y<- factor(y, levels=c(-1,1))
	y <- as.numeric(as.character(y))
	
	
	# scad svm model
	f<-model$f
	
	n = length(f)
	ind = rep(0, n)
	yf = y * f
	wt = ind + 2 * as.numeric(yf < -1) + as.numeric(yf <= 1 & yf >= -1)
	d = diag(model$xqx) * (0.5+1/abs(y-f))/n
	gacv = 0.5*mean((1-yf)+abs(1-yf))+mean(d*wt)
	return( gacv=gacv)
}

