wlai <- structure(function#Walter-Lieth aridity index
### Computing the annual aridity index from Walter-Lieth climate diagrams 
##details<<Areas between temperature and precipitation lines when precipitation exceeds temperature are calculated as indicators of moist seasons, and areas where temperature exceeds precipitation are calculated as indicator of dry season. The aridity index is defined as the quotient between the areas of dry and wet seasons. Those precipitations over 100 mm are scaled such that 1 degree C is equal to 5 mm.
##references<< Manrique E., A. Fernandez-Cancio. 2000. Extreme climatic events in dendroclimatic reconstructions from Spain. Clim. Chang., 44: 123-138.                  
(
    cd,##<< \code{data.frame} object with annual climatic records of
       ##monthly precipitation sums (mm), and monthly average
       ##temperatures (degree C), with row names being monthly
       ##characters in \code{\link{month.abb}} or
       ##\code{\link{month.name}}.
    sqt = TRUE, ##<<\code{logical}. Print the square root of the
                ##aridity index.  If TRUE then computed aridity index
                ##is normalized with a square root transformation.
    fig = FALSE, ##<<\code{logical}. Plot the Walter-Lieth diagram.
    ... ##<< Further arguments to be passed to \code{\link{plot}}
) {
    AI <- NA
    
    csn <- c(colclass(cd,T)[
                         c('tmp','fac')],recursive = T)
    csn. <- length(csn)!=0
    
    
    pr <- cd[,1];tm <- cd[,2]
    min.length <- length(pr) >= 12
    if(min.length){
        tm <- 2 * tm
        names(pr) <- 1:length(pr)
        names(tm) <- names(pr)
        
        fint <- function(pr,tm, twice = FALSE){
            if(twice)
                tm <- 2 * tm #twice temp
            ix<-which(diff(pr>tm)!=0)
            pr.m<-pr[ix+1]-pr[ix]
            tm.m<-tm[ix+1]-tm[ix]
            fx <- (tm[ix] - pr[ix]) / (pr.m-tm.m)
            ixc <- ix + fx
            prc <- pr[ix] + (pr.m*(ixc-ix))
            names(pr) <- names(tm) <- 1:length(pr)
            names(prc) <- ixc
            fx <- (tm[ix] - pr[ix]) / (pr.m-tm.m)
            ixc <- ix + fx
            prc <- pr[ix] + (pr.m*(ixc-ix))
            names(prc) <- ixc
            return(prc)}
        
        ford <- function(x){
            order(as.numeric(names(x)))}
        
        pr. <- ifelse(pr > 100, 80 + 0.2 * pr,pr)
        pr.. <- ifelse(pr > 100,100,pr)
        pics <- fint(rep(100,length(pr)),pr)#}
        prc <- fint(pr,tm)
        
        fxn <- function(x3,x4){
            xn <- c(ifelse(pr > tm,x3,x4),prc)
            ix <- order(as.numeric(names(xn)))
            xn <- xn[ix]
            return(xn)}
        
        xn <- fxn(pr.,tm)
        xn <- c(xn,pics)[ford(c(xn,pics))]
        xn1 <- fxn(tm,pr.)
        xn.. <- fxn(pr..,tm)
        xn.. <- c(xn..,pics)[ford(c(xn..,pics))]
        pr. <- c(pr.,pics)[ford(c(pr.,pics))]
        
        fb <- function(x){
            nx <- c(ptm,x,ptm)
            y <- as.numeric(names(x))
            names(nx) <- c(min(y),names(x),max(y))
            nxn <- as.numeric(names(nx)) - 0.5
            dt <- data.frame(x = nxn, y = nx)
            return(dt)}
        
        surf<-function(y){
            x <- y[,1]; y <- y[,2]
            dosarea<-sapply(2:(length(y)-1),
                            function(i)y[i]*(x[i+1]-x[i-1]))
            a <- 0.5 * sum(dosarea)
            return(a)}
        
        ptm <- min(pr.,tm)
        wet <- surf(fb(xn)) - surf(fb(tm))
        dry <- surf(fb(xn)) - surf(fb(pr.))
        AI <- dry/wet
        
        if(sqt)AI <- sqrt(AI)}
    
    if(csn.)
        AI <- unique(cbind(AI,cd[,csn[-1L]]))
    
    
    if(fig&min.length){
        ##color palette
        col. <- paste('gray',c(70,60,30), sep = '')
        
        par(oma = c(0,0,0,2))
        plot(fb(xn),
             ylim = c(ptm,max(c(pr.,tm))),
             col = col.[1],type='l',
             xaxt = 'n',yaxt = 'n',
             xlab = 'Month',ylab = NA, ...)
        mns <- rownames(cd)
        xap <- unique(fb(pr)[,1])
        yap <- pretty(c(ptm,max(c(pr.,tm))))
        axis(4,at = yap, labels = yap,las = 1)
        axis(2,at = yap, labels = yap/2,las = 1)
        axis(1,at = xap, tick = FALSE,labels = mns)            
        axis(1,at = xap - 0.5, tick = TRUE,labels = FALSE)            
        mtext(text = expression(~degree~C),
              las = 1,at = min(xap) - 1.2 * par('cex'))
        mtext(text = 'mm',las = 1,
              at = max(xap) + 1.2 *par('cex'))
        polygon(fb(xn), border = NA, col = col.[2])
        polygon(fb(xn..), border = NA, col = col.[1])
        polygon(fb(tm), border = NA, col = col.[3])
        polygon(fb(xn1), border = col.[1], col = 'white')
        lines(c(min(xap),max(xap)),c(100,100),col =col.[1],lty = 1)
        lp <- 'top'
        legend(lp,legend = c('Dry season','Moist season'),
               fill = c(col.[3],col.[1]),cex = 0.8,
               horiz = TRUE, border = NA,bty = 'n')}
    return(AI)                
### \code{numeric} aridity index and plot of the Walter-Lieth diagram.
} ,
ex=function()
{
    ##random simulation of climatic records
    set.seed(1)
    pr <- rnorm(12,1,1)
    tm <- rnorm(12,0,1)
    cld <- data.frame(pr,tm)
    ##labels of months from october to september
    rownames(cld) <- month.abb[c(10:12,1:9)]
    rownames(cld) <- c(10:12,1:9)
    ##computation of the aridity index and climate diagram
    AI <- wlai(cld,fig = TRUE)
    AI
})
