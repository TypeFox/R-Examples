plot.BIPOD <- function(x,
                       theta=NULL,
                       subset = NULL,
                       type,
                       lag = 20,
                       interval = NULL,
                       treshold = 0.1,
                       speed = 0.1,
                       truepath = NULL,
                       BY = 1,
                       prop=c(.05,.95),
                       diffPriorMean=NULL,
                       diffPriorCovar=NULL,
                       log=FALSE,
                       ...){
    if(!(type %in%
         c("trace","hist","acp","pairs","SDtrace","accept","movie","cover"))){
        stop("Invalid value of type. See help page for accepted values")}

    Ddrift1 <- dim(x$Drift)[1]
    Ddrift2 <- dim(x$Drift)[2]
    Ddiff1  <- dim(x$Diff)[1]
    Ddiff2  <- dim(x$Diff)[2]
    Ndrift  <- ShowModels(x$Info$Model)$Ndrift
    Ndiff   <-ShowModels(x$Info$Model)$Ndiff

    Names2 <- c(paste("drift",1:Ndrift,sep=""),
                paste("diff", 1:Ndiff,sep=""))

    if(is.null(interval)){interval <- 1:max(Ddrift1,Ddiff1,dim(x$LatentPath)[2])}
    if(length(interval)==1){
        interval <- interval:max(Ddrift1,Ddiff1,dim(x$LatentPath)[2])
    }
    if(is.null(subset)){
        subset <- 1:(Ndrift+Ndiff)
        Nam <- Names2[subset]
    }

    if(Ddrift1==1 & Ddiff1==1 & !(type %in% c("movie","cover","accept"))){
        stop("This BIPOD object has no drift or diffusion estimates to plot")
    }
    if(Ddrift1!=1 & Ddiff1==1){
      subset <- intersect(subset,1:Ndrift)
      Nam <- Names2[subset]
      xx <- x$Drift[interval,subset]
      if(length(subset)==0){
        stop(paste("Argument 'subset' must be a subset of 1:",Ndrift,sep=""))
      }
    }
    if(Ddrift1==1 & Ddiff1!=1){
        subset <- intersect(subset-Ndrift,1:Ndiff)
        Nam <- Names2[subset+Ndrift]
        xx <- x$Diff[interval,subset]
        if(length(subset)==0){
          stop(paste("Argument 'subset' must be a subset of ",Ndrift+1,":",Ndrift+Ndiff,sep=""))
        }
      }
    if(Ddrift1!=1 & Ddiff1!=1){
        Nam <- Names2[subset]
        xx <- cbind(x$Drift,x$Diff)[interval,subset]
    }
    if(length(subset)==1){xx <- as.matrix(xx)}


    par(mfrow = c(2, floor((length(subset) + 1)/2)),las=1)
    if (length(subset) == 1) {
        par(mfrow = c(1, 1),las=1)
    }


    if (type == "trace") {
        for (i in 1:length(subset)) {
            matplot(interval,
                    xx[,i],
                    type = "l",
                    col = 1,
                    xlab="",
                    ylab="",
                    main = Nam[i],
                    ...)
            if(!is.null(theta)){
                abline(h = theta[subset[i]], col = "gray")
            }
        }
    }

    if (type == "hist") {
        PriorSpec <- TRUE
        if(is.null(diffPriorMean) & is.null(diffPriorCovar)){
            PriorSpec <- FALSE}
        for (i in 1:length(subset)) {
            if(subset[i]==5 & log==TRUE){
                xx[,i] <- log(xx[,i])
                theta[5] <- log(theta[5])
            }
            if(subset[i]==6 & log==TRUE){
                xx[,i] <- log(xx[,i])
                theta[6] <- log(theta[6])
            }

            plot(density(xx[, i],...),
                 main = Nam[i],
                 xlab="")
            if(!is.null(theta)){
                abline(v = theta[subset[i]],
                       col = "gray",
                       lwd = 3)
            }
            if(subset[i] == 5){
                if(log==FALSE & PriorSpec==TRUE){
                    curve(dlnorm(x,meanlog=diffPriorMean[1],sdlog=sqrt(diffPriorCovar[1,1])),
                          add=T,
                          col="gray",
                          lty=2)
                }
                if(log==TRUE  & PriorSpec==TRUE)
                    curve(dnorm(x,mean=diffPriorMean[1],sd=sqrt(diffPriorCovar[1,1])),
                          add=T,
                          col="gray",
                          lty=2)
            }
        }

      if(subset[i] == 6 & PriorSpec==TRUE){
          curve(dnorm(x,mean=diffPriorMean[2],sd=sqrt(diffPriorCovar[2,2])),
              add=T,
              col="gray",
              lty=2)
    }
  }

        if (type == "acp") {
        for (i in 1:length(subset)) {
            acf(xx[, i],
                main = Nam[i],
                lag.max = lag,
                ...)
        }
    }
    if (type == "pairs") {
      pairs(xx,...)
    }
    if (type == "SDtrace") {
        standfun <- function(x){
            (x - mean(x))/sd(x)
        }
        xx2 <- apply(xx, 2, standfun)
        for (i in 1:length(subset)) {
            matplot(xx2[, i],
                    type = "l",
                    col = 1,
                    ylab = "",
                    main = Nam[i],
                    ...)
            abline(h = 0,
                   col = 1)
            if(!is.null(theta)){
                abline(h = (theta[i] - mean(xx[, i]))/sd(xx[,i]),
                       col = "gray")
            }
        }
    }
    if (type == "accept") {
        par(mfcol = c(2, 2), las = 1)
        plot(yy <- apply(x$PathAcc, 2, mean),
             type = "l",
             xlab = "Iterations",
             ylim = range(c(yy, treshold)),
             ylab = "Acc rate",
             main="Path Acc")
            abline(h = 0:10/10, col = "lightgray", lty = "dotted")

        dd <- which(apply(x$PathAcc, 1, mean) < treshold)
        plot(apply(x$PathAcc, 1, mean),
             type = "l",
             main = paste("Path;", length(dd),
               " timepoints with acc. rate < ",
               treshold, sep = ""),
             xlab = "Latent data paths",
             ylab = "Acc rate")
        points(dd, rep(treshold, length(dd)), pch = 19, col = 2)

        plot(yy <- apply(x$PathAccPoints, 2, mean),
             type = "l",
             xlab = "Iterations",
             ylim = range(c(yy, treshold)),
             ylab = "Acc rate",
             main="Point Acc")
        abline(h = 0:10/10, col = "lightgray", lty = "dotted")

        dd <- which(apply(x$PathAccPoints, 1, mean) < treshold)
        plot(apply(x$PathAccPoints, 1, mean),
             type = "l",
             main = paste("Point;", length(dd),
               " timepoints with acc. rate < ",
               treshold, sep = ""),
             xlab = "Latent data points",
             ylab = "Acc rate")
        points(dd, rep(treshold, length(dd)), pch = 19, col = 2)

      }

    if (type == "movie"){
        if(x$Info$Coords=="Two coordinates observed"){
          stop("Plot type not valid when both coordinates are observed\n")}
        par(mfrow=c(1,1),las=1)
        ## for(i in seq(1,dim(x$LatentPath)[2],by=BY)){
        for(i in seq(interval[1],tail(interval,n=1),by=BY)){
            ##for(i in interval){
            matplot(x$LatentPath[,i],
                    type="l",
                    col=1,
                    ylim=range(x$LatentPath,truepath),
                    main=i,
                    xlab="Observed timepoints",
                    ylab="Latent path values")
            if(!is.null(truepath)){matlines(truepath,col=2,lwd=2,type="l")}
            Sys.sleep(speed)
        }
    }

    if (type == "cover"){
        if(x$Info$Coords=="Two coordinates observed"){
            stop("Plot type not valid when both coordinates are observed\n")
        }
        par(mfrow=c(1,1),las=1)
        tmp <- apply(x$LatentPath[,interval],1,quantile,p=c(prop,.5))
        tmp2 <- apply(x$LatentPath[,interval],1,mean)
        yy <- range(tmp,truepath)
        matplot(NA,
                xlim=c(1,length(truepath)),
                ylim=yy,
                ylab="",
                las=1,
                ...)
#                main=paste(prop[1],"-",prop[2],"quantiles",sep=""))
        for(i in 1:(dim(x$LatentPath)[1])){
         lines(rep(i,2),tmp[1:2,i],col="gray",lwd=2)
#         points(rep(i,2),c(tmp[3,i],tmp2[i]),pch=19,col=1:2)
       }
       matlines(truepath,type="l",...)
     }
}
