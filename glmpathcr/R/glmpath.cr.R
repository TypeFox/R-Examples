glmpath.cr <-
function(x,y,data,method="backward",weight = rep(1, n), offset = rep(0, n), lambda2 = 1e-05, 
    max.steps = 10 * min(n, m), max.norm = 100 * m, min.lambda = (if (m >= 
        n) 1e-06 else 0), max.vars = Inf, max.arclength = Inf, 
    frac.arclength = 1, add.newvars = 1, bshoot.threshold = 0.1, 
    relax.lambda = 1e-08, standardize = TRUE, function.precision = 3e-13, 
    eps = .Machine$double.eps, trace = FALSE) {
	n<-dim(x)[1]
	p<-m<-dim(x)[2]
	k<-length(unique(y))
	x<-as.matrix(x)
	if (c("backward","forward")[charmatch(method,c("backward","forward"))]=="backward") {
		restructure<-cr.backward(x=x,y=y)
	}
	if (c("backward","forward")[charmatch(method,c("backward","forward"))]=="forward") {
		restructure<-cr.forward(x=x,y=y)
	}		
   glmpath.data<-list(x=restructure[,-1],y=restructure[,"y"])
   object<-glmpath(x, y, data=glmpath.data, family=binomial, standardize=TRUE, nopenalty.subset=(p+1):(p+k-1))
   object$x<-x
   object$y<-y
   object$method<-method
   class(object)<-"glmpath.cr"
   object
}

