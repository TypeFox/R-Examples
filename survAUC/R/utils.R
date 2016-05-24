


my.sort <- function(x,index.return=FALSE,decreasing=FALSE) {
	if(decreasing) 
		x <- -x
	y <- sort(x, method="quick", index.return=index.return)
	if(decreasing) 
		y$x <- -y$x
	y
}


################################################################
###				step function                    			 ###
################################################################
## x		- x-Werte
## y		- y-Werte
## z		- neue x-Werte



seval <- function(x,y,z){
	n_x <- length(x)
	n_z <- length(z)
	xz_new <- rep(0,n_x) 
	ans <- .C("step_eval_R",
			  as.numeric(xz_new),
			  as.numeric(z), 
			  as.numeric(x), 
			  as.numeric(y), 
			  as.integer(n_z),
			  as.integer(n_x),
			  PACKAGE="survAUC")
	ans[[1]]
}