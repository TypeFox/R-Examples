amod <- structure(function#Allometric modeling
### Parameters in simple allometric model are evaluated on tree radial increments to compute diameters, basal areas, tree biomasses, etc.
##details<<The simple allometric model has the form: a * cs ^ b, with a,b being constants in \code{mp}, and cs being scaled-cummulative sums.  Different dendrometric variables can be computed; for example, \code{c(1,1)} produces diameters, and \code{c(0.25 * pi,2)} computes basal areas. Argument \code{mp} can have more than two parameters: \code{c(a1,b1,a2,b2, ..., an,bn)}, with \code{n} being the number of times that allometric model will be recursively implemented. Such recursive evaluation is useful to derive variables which depend on other allometric covariables: i.e allometric model would be implemented twice to recursively compute diameters and tree biomasses. A column of increments of cs (x) is also computed by implementing \code{\link{setdiff}}. 
(
    
    cs, ##<<\code{Numeric} vector of scaled-cummulative sums such as
    ##that produced by \code{\link{scacum}}.
    mp = c(0.5,1), ##<<\code{Numeric}. vector with allometric
    ##parameters. Default \code{c(0.5,1)} maintains the
    ##original radii (see details for other variables)
    un = NULL ##<< NULL, or bidimensional \code{character} vector to
    ##transform SI units of the processed variable. The SI
    ##units can be expressed in micrometers 'mmm',
    ##milimeters 'mm', centimeters 'cm', decimeters 'dm', or
    ##meters 'm'. If NULL then original units are
    ##maintained.
) {
    
    csn. <- FALSE
    if(is.data.frame(cs)){
        csnu <- colclass(cs,T)[['num']]
        csn <- c(colclass(cs,T)[
                             c('tmp','fac')],recursive = T)
        csn. <- length(csn)!=0
        csn.. <- csn[!csn%in%c('x','csx')]
        cd <- cs
        cs <- cs[,'csx']
        names(cs) <- cd[,'year']}
    
    chun <- function(from,to){
        sm <- 10 ^ -c(6,3:0)
        un <- c('mmm','mm','cm','dm','m')
        names(sm) <- un
        eq <- sm[from]/sm[to]
        names(eq) <- to
        return(eq)}
    allm <- function(x,a,b){
        a * (x ^ b)}
    x <- 2 * cs #diameters
    if(length(un) > 1)
        x <- x * chun(un[1],un[2])
    
    fp <- function(x){
        cn <- c(TRUE,FALSE)
        l <- list(a = x[cn],b = x[rev(cn)])
        return(l)}
    
    x0 <- x
    arg <- list()
    for(i in 1:length(fp(mp)[['a']])){
        arg[[i]] <- c(fp(mp)[['a']][i],
                      fp(mp)[['b']][i])
        x0 <- do.call(allm,list(x0,
                                arg[[i]][1],arg[[i]][2]))}
    
    x1 <- c(NA,diff(x0))
    names(x1) <- names(x)
    xd <- data.frame(x = x1,csx = x0)
    
    if(csn.&& length(csnu) > 1){
        xd <- cd[,csnu]
        xd[,'x'] <- x1
        xd[,'csx'] <- x0 }
    if(csn.)
        xd <- cbind(xd,cd[,csn..])
    
    return(xd)
### \code{data.frame} object.
    
} ,
ex=function() {
    ## radial increments
    set.seed(1)
    w <- abs(rnorm(12,1,1))
    names(w) <- 1951:1962
    ## scaled and cummulative radial increments
    sr <- scacum(w)
    ## diameters
    d <- amod(sr[,2],c(1,1))
    ## basal areas (m2):
    ba <- amod(sr[,2],c(0.25 * pi,2),c('mm','m'))
    print(ba)
})
