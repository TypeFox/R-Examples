matplotLong <- function (trajLong, col = 1:6,lty=1:5,lwd = 1, add = FALSE, main="",xlab="Times",ylab="", pourcent=NA){
    if (!add) {
        xlimit <- range(trajLong[, 2], na.rm = TRUE)
        ylimit <- range(trajLong[, 3], na.rm = TRUE)
        plot(1, type = "n", ylim = ylimit, xlim = xlimit, main=main,xlab=xlab,ylab=ylab)

    }else{}

    if(!identical(pourcent,NA)){
        legend("top",legend=paste0(round(pourcent*100),"%"),lty=1,pch="",
               ncol=legendCol(length(pourcent)),xpd=NA,xjust=0.5,inset=-0.12,col=1:length(pourcent)+1)
    }else{}

    indiv <- factor(unique(trajLong[, 1]))
    if (length(col) < length(indiv)) {
        col <- rep_len(col, length(indiv))
    }else{}
    names(col) <- indiv

    if (length(lty) < length(indiv)) {
        lty <- rep_len(lty, length(indiv))
    }else{}
    names(lty) <- indiv

    for (iIndiv in indiv) {
        select <- trajLong[, 1] == iIndiv
        lines(x = trajLong[select, 2], y = trajLong[select, 3],col = col[iIndiv], lwd = lwd,lty=lty[iIndiv])
    }
    return(invisible())
}

matplotWide <- function (trajWide, times=1:ncol(trajWide), col = 1:6,lty=1:5,lwd = 1, main="",xlab="Times",ylab="", pourcent=NA){
    matplot(times,t(trajWide),col=col,lty=lty,lwd=lwd,main=main,xlab=xlab,ylab=ylab,type="l")

    if(!identical(pourcent,NA)){
        legend("top",legend=paste0(round(pourcent*100),"%"),lty=1,pch="",
               ncol=legendCol(length(pourcent)),xpd=NA,xjust=0.5,inset=-0.12,col=1:length(pourcent)+1)
    }else{}

    return(invisible())
}


## matplotLongLegend <- function (trajLong, col = 1:6, lty=1:5, lwd = 1, pourcent,...){
##     xlimit <- range(trajLong[, 2], na.rm = TRUE)
##     ylimit <- range(trajLong[, 3], na.rm = TRUE)
##     plot(1, type = "n", ylim = ylimit, xlim = xlimit,...)
##     legend(x=mean(xlimit),y=ylimit[2]*1.21-ylimit[1]*0.21,legend=paste0(round(pourcent*100),"%"),lty=1,pch="",
##            ncol=legendCol(length(pourcent)),xpd=NA,xjust=0.5,inset=-0.1,col=1:length(pourcent)+1)

##     indiv <- factor(unique(trajLong[, 1]))
##     if (length(col) < length(indiv)) {
##         col <- rep_len(col, length(indiv))
##     }else{}
##     names(col) <- indiv

##     if (length(lty) < length(indiv)) {
##         lty <- rep_len(lty, length(indiv))
##     }else{}
##     names(lty) <- indiv

##     for (iIndiv in indiv) {
##         select <- trajLong[, 1] == iIndiv
##         lines(x = trajLong[select, 2], y = trajLong[select, 3],col = col[iIndiv], lwd = lwd)
##     }
##     return(invisible())
## }



ClusterLongDataShape_plotTraj <- function(x,y,col="clusters",pourcent=NA,...){
    if(identical(col,"clusters")){
        if(x["kmlShape"]){
            col <- as.integer(x["clusters"])+1
        }else{
            col <- "grey"
        }
    }else{}

    if(x["longAvailable"]){
       matplotLong(x["trajLong"],col=col,pourcent=pourcent,...)
    }else{
       matplotWide(trajWide=x["trajWide"],times=x["times"],col=col,pourcent=pourcent,...)
    }
    return(invisible())
}
setMethod("plotTraj",signature=c("Clds","missing"),ClusterLongDataShape_plotTraj)



ClusterLongDataShape_plotMeans <- function(x,y,add=FALSE,pourcent=NA){
    matplotLong(x["trajMeans"],col=1, lwd=8,lty=1,add=add,pourcent=pourcent)
    matplotLong(x["trajMeans"],col=2:(x["nbClusters"]+1), lwd=4,lty=1,add=TRUE,pourcent=pourcent)
    return(invisible())
}
setMethod("plotMeans",signature=c("Clds","missing"),ClusterLongDataShape_plotMeans)


ClusterLongDataShape_plot <- function(x,y,col="darkgrey",lty=1,legend=TRUE,...){
    if(legend & x["kmlShape"]){
        plotTraj(x,col=col,lty=lty,pourcent=table(x["clusters"])/length(x["clusters"]),...)
    }else{
        plotTraj(x,col=col,lty=lty,...)
    }
    if(x["kmlShape"]){
            plotMeans(x,add=TRUE)
    }else{}
    return(invisible())
}
setMethod("plot",signature=c("Clds","missing"),ClusterLongDataShape_plot)




plotSenators <- function(x,col = 2:7,lty=1:5,lwd = 1, add = FALSE, main="",xlab="Times",ylab=""){
    if(x["senatorsAvailable"]){
        matplotLong(x["senators"],col=col,lty=lty,lwd=lwd,add=add,main=main,xlab=xlab,ylab=ylab)
    }else{
        warning("[kmlShape:plotSenator] There is no senator to plot !")
    }
    return(invisible())
}




