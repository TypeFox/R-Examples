`accumcomp` <-
function(x,y="",factor,scale="",method="exact",permutations=100,conditioned=T,gamma="boot",plotit=T,labelit=T,legend=T,rainbow=T,xlim=c(1,max),ylim=c(0,rich),type="p",xlab="sites",ylab="species richness",...) {
    groups <- table(y[,factor])
    min <- min(groups)
    max <- max(groups)
    m <- length(groups)
    levels <- names(groups)
    result <- array(NA,dim=c(m,max,3))
    dimnames(result) <- list(level=levels,obs=c(1:max),c("Sites","Richness","sd"))
    names(dimnames(result)) <- c(factor,"obs","")
    for (i in 1:m) {
        result1 <- accumresult(x,y,factor,level=levels[i],scale=scale,method=method,permutations=permutations,conditioned=conditioned,gamma=gamma)
        l <- length(result1$sites)
        result[i,c(1:l),1] <- result1$sites
        result[i,c(1:l),2] <- result1$richness
        if (method!="collector" && method!="poisson" && method!="binomial" && method!="negbinomial") {result[i,c(1:l),3] <- result1$sd}
    }
    if (plotit == T) {
        max <- max(result[,,1],na.rm=T)
        rich <- max(result[,,2],na.rm=T)
        for (i in 1:m) {
            result1 <- accumresult(x,y,factor,level=levels[i],scale=scale,method=method,permutations=permutations,conditioned=conditioned,gamma=gamma)
            if (plotit == T) {
                if (i == 1) {addit <- F}
                if (i > 1) {addit <- T}
                if (labelit==T) {
                    labels <- levels[i]
                }else{
                    labels <- ""
                }
                if (rainbow==T) {
                    grDevices::palette(rainbow(m))
                    accumplot(result1,addit=addit,xlab=xlab,ylab=ylab,xlim=xlim,ylim=ylim,labels=labels,col=i,pch=i,type=type,...)
                }else{
                    accumplot(result1,addit=addit,xlab=xlab,ylab=ylab,xlim=xlim,ylim=ylim,labels=labels,pch=i,type=type,...)
                }
            }
        }
        if (rainbow==T && legend==T) {legend(graphics::locator(1),legend=levels,pch=c(1:m),col=c(1:m))}
        if (rainbow==F && legend==T) {legend(graphics::locator(1),legend=levels,pch=c(1:m))}
    }
    grDevices::palette("default")
    return(result)
}

