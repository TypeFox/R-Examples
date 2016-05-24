Xgen <- function(tt,n,s, pp, design,  steps=10) {
    # Builds design matrix
    # Arguments not checked as it is assumed that the function singlePer gives the right arguments
    # Fold time points to phase [0,1]
    t_ <- (tt%%pp)/pp
    
    # minimum points that are needed to fit a function
    if(design%in%c("step", "stepB")) points.needed <- 1
    if(design=="sine") points.needed <- 3
    if(design=="fourier(2)") points.needed <- 5
    if(design=="fourier(3)") points.needed <- 7
    if(design=="splines") points.needed <- 4
    
    if(length(unique(t_))<points.needed) {
        warning(paste("trial period", pp, "could not be fitted due to improper sampling"))
        return(NA)
    }
    
    # Build designmatrix
    if(design=="step") {
        level <- 1+(t_*steps)%/%1
        X <- matrix(0, nrow=n, ncol=steps)
        X[cbind(1:n,level)] <- 1 
        X <- X[,which(apply(X,2,max)!=0)]
    }
        
    if(design=="stepB") {
        level <- 1+((0.5+t_*steps)%/%1%%steps)
        X <- matrix(0, nrow=n, ncol=steps)
        X[cbind(1:n,level)] <- 1 
        X <- X[,which(apply(X,2,max)!=0)]
    }
    
    if(design%in%c("sine", "fourier(2)", "fourier(3)")) {
        X <- cbind(1,sin(t_*2*pi), cos(t_*2*pi))
    }
    
    if(design%in%c("fourier(2)", "fourier(3)")) {
        X<- cbind(X, sin(t_*4*pi), cos(t_*4*pi))
    }
    
    if(design=="fourier(3)") {
        X <- cbind(X, sin(t_*6*pi), cos(t_*6*pi))
    }
       
    if(design=="splines") {
        temp<-rbind(diag(4), cbind(diag(3), 0))
        X <- spline.des(knots=seq(from=-0.75, to=1.75, by=0.25), t_)$design%*%temp
    }   
	
    # Weight and return design matrix
    X.weighted <- cbind(X/s)
    
    if(qr(X.weighted)$rank<points.needed) {
        warning(paste("trial period", pp, "could not be fitted due to improper sampling"))
        return(NA)
    } else return(X.weighted)
}
