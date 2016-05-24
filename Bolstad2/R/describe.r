describe<-function(x, varNames = NULL){
    ## Mimics Minitab's desc function

    nameX<-deparse(substitute(x))

    if(is.matrix(x)){
        x<-data.frame(x)
        if(is.null(varNames))
            varNames<-paste(nameX,0:(ncol(x)-1),sep="")
        names(x)<-varNames
    }

    nx<-sapply(x,length)
    mx<-sapply(x,mean)
    sx<-sapply(x,sd)
    SEx<-sapply(x,sd)/sqrt(nx)
    minx<-sapply(x,min)
    maxx<-sapply(x,max)
    q1x<- sapply(x,quantile,prob=0.25)
    medx<-sapply(x,median)
    q3x<-sapply(x,quantile,prob=0.75)

    stats.df<-data.frame(N=nx, mean = mx, stdev = sx, sterr = SEx,
                        min = minx, q1 = q1x, med = medx, q3 = q3x, max = maxx)

    nVars <- ncol(x)

    print(stats.df)

    invisible(stats.df)
}
