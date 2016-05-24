`rankabuncomp` <-
function(x,y="",factor,scale="abundance",scaledx=F,type="o",rainbow=T,legend=T,xlim=c(1,max1), ylim=c(0,max2), ...) {
    groups <- table(y[,factor])
    levels <- names(groups)
    m <- length(groups)
    max1 <- max(diversitycomp(x,y,factor,index="richness")[,2])
    if (scaledx==T) {xlim<-c(0,100)}
    freq <- diversityresult(x,index="Berger")
    if (scale=="abundance") {max2 <- freq * diversityresult(x,index="abundance")}
    if (scale=="logabun") {max2 <- log(freq * diversityresult(x,index="abundance"),base=10)}
    if (scale=="proportion") {max2 <- 100 * max(diversitycomp(x,y,factor,index="Berger")[,2])}
    if (scale=="accumfreq") {max2 <- 100}
    max2 <- as.numeric(max2)
    if (rainbow==F) {
        rankabunplot(rankabundance(x, y, factor, levels[1]), scale=scale, scaledx=scaledx, type=type, labels=levels[1], xlim=xlim, ylim=ylim, pch=1, specnames=NULL, ...)
        for (i in 2:m) {
            rankabunplot(rankabundance(x, y, factor, levels[i]), addit=T, scale=scale, scaledx=scaledx, type=type, labels=levels[i], pch=i, specnames=NULL,...)
        }
        if (legend==T) {legend(graphics::locator(1),legend=levels,pch=c(1:m))}
    }else{
        grDevices::palette(rainbow(m))
        rankabunplot(rankabundance(x, y, factor, levels[1]), scale=scale, scaledx=scaledx, type=type, labels=levels[1], xlim=xlim, ylim=ylim, col=1, pch=1, specnames=NULL,...)
        for (i in 2:m) {
            rankabunplot(rankabundance(x, y, factor, levels[i]), addit=T, scale=scale, scaledx=scaledx, type=type, labels=levels[i], col=i, pch=i, specnames=NULL,...)
        }
        if (legend==T) {legend(graphics::locator(1), legend=levels, pch=c(1:m), col=c(1:m))}
        grDevices::palette("default")
    }
}

