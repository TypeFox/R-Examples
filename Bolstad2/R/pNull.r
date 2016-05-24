pNull<-function(theta0, theta, cdf = NULL, type = 'upper'){

    if(length(theta)<10)
        stop("theta must have at least ten values")

    if(length(grep('^[lL]',type))>0){
        type<-'lower'
    }else if(length(grep('^[Uu]',type))>0){
        type<-'upper'
    }else{
        stop("type must be one of lower or upper")
    }

    Fx<-ecdf(theta)

    if(!is.null(cdf))
    {
        o <- order(theta)
        if(any(theta[o]!=theta)){
            warning("theta is not in ascending order. This may cause problems")
        }
        suppressWarnings(Fx<-approxfun(theta,cdf))
    }

    if(type=='lower'){
        prob<-Fx(theta0)
        cat(paste("Posterior Pr(theta<=theta0) is ",
                  prob, "\n",sep=""))
        invisible(list(prob=prob))
    }else{
        prob<-1-Fx(theta0)
        cat(paste("Posterior Pr(theta>=theta0) is ",
                  prob, "\n",sep=""))
        invisible(list(prob=prob))
    }
}
