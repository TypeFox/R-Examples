pnullNum<-function(theta0, theta, cdf, type = 'upper'){

    if(length(theta)<10)
        stop("theta must have at least ten values")

    if(length(grep('^[lL]',type))>0){
        type<-'lower'
    }else if(length(grep('^[Uu]',type))>0){
        type<-'upper'
    }else{
        stop("type must be one of lower or upper")
    }

    Fx<-approxfun(theta,cdf)

    if(type=='lower'){
        prob<-1-Fx(theta0)
        cat(paste("Posterior Pr(theta<=theta0) is ",
                  prob, "\n",sep=""))
        invisible(list(prob=prob))
    }else{
        prob<-Fx(theta0)
        cat(paste("Posterior Pr(theta>=theta0) is ",
                  prob, "\n",sep=""))
        invisible(list(prob=prob))
    }
}
