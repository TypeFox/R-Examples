GANOVA <- function(dataset, var.equal=TRUE, type="QQ"){
    dataset <- dataset[order(dataset[,1]),]
    data.split <- split(dataset[,2], dataset[,1])
    treatments <- sapply(data.split,mean)
    n <- sapply(data.split, length)
    errors <- dataset[,2]-rep(treatments, n)
    if (var.equal) {
        treatments <- (treatments-mean(treatments))*sqrt(n)
    } else {
        sigmas <- sapply(data.split, sd)
        errors <- errors/rep(sigmas, n)
        treatments <- (treatments-mean(treatments))*sqrt(n)/sigmas
    }
    trt.as.factor <- factor(dataset[,1])
    plotregionmin <- min(treatments, errors)
    plotregionmax <- max(treatments, errors)
    if (type=="QQ") {
        qqANOVA(errors, treatments, xlim=range(c(treatments, errors)), 
           ylim=range(c(treatments, errors)), ylab="scaled treatment averages")
        abline(0,1)
        legend(plotregionmin, plotregionmax, legend=levels(trt.as.factor), col=as.numeric(factor(treatments)), pch=16)
        tickhalflength <- sd(treatments)/10
        segments(errors, rep(plotregionmin, length(errors))-tickhalflength, errors, rep(plotregionmin, length(errors))+tickhalflength, 
            col=as.numeric(trt.as.factor), lwd=2)
    } else {    
        breaks <- hist(c(treatments, errors), breaks="FD", plot=FALSE)$breaks
        hist(errors, xlab="errors", main="", breaks=breaks)   
        points(treatments, rep(0, length(treatments)), cex=2.5, col=1:length(treatments),
           pch=7+1:length(treatments))
        legend("topleft", legend= c(paste("average", unique(dataset[,1]))),
           pch=7+1:length(treatments), col=1:length(treatments))
    }
}
