plot.PredictFindIt <- function(x,main,xlab,ylab, ...){
    labels <- "index"
    if(missing(xlab)){
        xlab <- "index of observation"
    }
    if(missing(ylab)){
        ylab <- "Treatment Effect"
    }
    object <- x
    treat.type <- object$treat.type
    if(treat.type=="single"){
        pred.data <- object$data
        ATE <- object$ATE
        pred.data.out.p <- pred.data[order(pred.data$Treatment.effect,
                                           decreasing=FALSE),]   
        xp <- seq(1:nrow(pred.data.out.p))
        zero.1 <- min(which(min(abs(pred.data.out.p$Treatment.effect))==
                            abs(pred.data.out.p$Treatment.effect)))
        zero <- xp[zero.1]
        low <- min(pred.data.out.p$Treatment.effect)
        if(missing(main)){
            main <- "Causal Moderation: Heterogeneous Treatment Effect "
        }
        plot(xp,pred.data.out.p$Treatment.effect,type="l",col="red",
             main = main,
             ylab= ylab,xlab=xlab)
        text(zero,low, labels=as.character(zero))
        abline(h=0,lty="dotdash")
        abline(h=ATE, col="blue")
        abline(v=zero, col="grey")                                    
        if(labels=="index"){
            p <- try(identify(xp,pred.data.out.p$Treatment.effect,
                              labels=rownames(pred.data.out.p)
                              ),silent=TRUE)
        }
        else{
            p <- try(identify(xp,pred.data.out.p$Treatment.effect,
                              labels=pred.data.out.p[,labels]
                              ),silent=TRUE)
        }
    }
    if(treat.type=="multiple"){
        pred.data <- object$data
        ## pred.data <- data[,3:ncol(data)]
        ATE <- object$ATE
        pred.data.out.p <- pred.data[order(pred.data$Treatment.effect,
                                           decreasing=FALSE),]   
        xp <- seq(1:nrow(pred.data.out.p))
        zero.1 <- min(which(min(abs(pred.data.out.p$Treatment.effect))==
                            abs(pred.data.out.p$Treatment.effect)))
        zero <- xp[zero.1]
        low <- min(pred.data.out.p$Treatment.effect)
        if(missing(main)){
            main <- "Causal Interaction: Heterogeneous Treatment Effect"
        }
        plot(xp,pred.data.out.p$Treatment.effect,col="red",
             main = main,
             ylab= ylab,xlab=xlab)
        text(zero,low, labels=as.character(zero))
        abline(h=0,lty="dotdash")
        abline(h=ATE, col="blue")
        abline(v=zero, col="grey")
        
        if(labels=="index"){
            p <- try(identify(xp,pred.data.out.p$Treatment.effect,
                              labels=rownames(pred.data.out.p)
                              ),silent=TRUE)
        }
        else{
            p <- try(identify(xp,pred.data.out.p$Treatment.effect,
                              labels=pred.data.out.p[,labels]
                              ),silent=TRUE)
        }
    }
    
}
