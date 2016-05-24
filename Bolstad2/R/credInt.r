credInt<-function(theta, cdf = NULL, conf = 0.95, type="twosided"){

    if(conf <=0  | conf >=1)
        stop("conf must be between 0 and 1")

    if(length(grep('^[Ll]',type))>0){
        type<-'lower'
    }else if(length(grep('^[Tt]',type))>0){
        type<-'twosided'
    }else if(length(grep('^[Uu]',type))>0){
        type<-'upper'
    }else{
        stop("type must be one of lower, upper or twosided")
    }

    alpha<-1-conf
    n<-length(theta)

    if(n<10)
        stop("theta must have at least ten values")

    if(!is.null(cdf)){
        suppressWarnings(Finv<-approxfun(cdf, theta))

        if(type=='lower'){
            lower.bound = Finv(alpha)
            cat(paste("Lower credible bound is : ", lower.bound, "\n",sep=""))
            invisible(list(lower.bound=lower.bound))
        }else if(type=='upper'){
            upper.bound<-Finv(1-alpha)
            cat(paste("Upper credible bound is : ", upper.bound, "\n",sep=""))
            invisible(list(upper.bound=upper.bound))
        }else{
            lower.bound = Finv(alpha/2)
            upper.bound<-Finv(1-alpha/2)
            cat(paste("Credible interval is : (", lower.bound,
                      ",", upper.bound,")\n",sep=""))
            invisible(list(lower.bound=lower.bound,upper.bound=upper.bound))
        }
    }else{
        if(type=='lower'){
            lower.bound <- quantile(theta, alpha)
            cat(paste("Lower credible bound is ",
                      lower.bound, "\n", sep=""))
            invisible(list(lower.bound=lower.bound))
        }else if(type=='upper'){
            upper.bound <- quantile(theta, 1-alpha)
            cat(paste("Upper credible bound is ",
                      upper.bound, "\n", sep=""))
            invisible(list(upper.bound=upper.bound))
        }else{
            bounds<-quantile(theta, c(alpha/2,1-alpha/2))
            lower.bound <-bounds[1]
            upper.bound <-bounds[2]
            cat(paste("Credible interval is (",
                      lower.bound, ',', upper.bound,")\n", sep=""))
            invisible(list(lower.bound=lower.bound,upper.bound=upper.bound))
        }
    }
}





