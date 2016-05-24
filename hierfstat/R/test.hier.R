"test.between" <-
function (data = data, test.lev, rand.unit, nperm = 100, ...) 
{
    get.g <- function(x, data, ...) {
        g.stats.glob(data.frame(x, data), ...)$g.stats
    }
    test.lev<-as.integer(factor(test.lev))
    rand.unit<-as.integer(factor(paste(test.lev,"X",rand.unit,sep="")))
    x<-order(test.lev,rand.unit)
    data<-data.frame(data[x,])
    test.lev<-test.lev[x]
    rand.unit<-rand.unit[x]
#    runit<-paste(test.lev,rand.unit,sep="")
#    runit<-rep(1:dim(table(runit)),c(table(runit)))
#    nobs <- length(runit)
    perm.stat <- vector(length = nperm)
    perm.stat[nperm] <- get.g(test.lev, data, ...)
    for (i in 1:(nperm - 1)) {
        perm.stat[i] <- get.g(test.lev, data[samp.between(rand.unit), 
            ], ...)
    }
    list(g.star = perm.stat, p.val = sum(perm.stat >= perm.stat[nperm])/nperm)
}
"test.between.within" <-
function (data = data, within, test.lev, rand.unit, nperm = 100, 
    ...) 
{
    get.g <- function(x, data, ...) {
        g.stats.glob(data.frame(x, data), ...)$g.stats
    }
    within<-as.integer(factor(within))
    test.lev<-as.integer(factor(paste(within,"X",test.lev,sep="")))
    rand.unit<-as.integer(factor(paste(within,"X",test.lev,"Y",rand.unit,sep="")))
    x<-order(within,test.lev,rand.unit)
    data<-data.frame(data[x,])
    within<-within[x]
    test.lev<-test.lev[x]
    rand.unit<-rand.unit[x]
#    tlev<-paste(within,test.lev,sep="")#
#    runit<-paste(tlev,rand.unit,sep="")#
#    tlev<-rep(1:dim(table(tlev)),c(table(tlev)))
#    runit<-rep(1:dim(table(runit)),c(table(runit)))
#    nobs <- length(test.lev)
    perm.stat <- vector(length = nperm)
    perm.stat[nperm] <- get.g(test.lev, data, ...)
    for (i in 1:(nperm - 1)) {
        perm.stat[i] <- get.g(test.lev, data[samp.between.within(inner.lev = rand.unit, 
            outer.lev = within), ], ...)
    }
    list(g.star = perm.stat, p.val = sum(perm.stat >= perm.stat[nperm])/nperm)
}
"test.g" <-
function (data = data, level, nperm = 100, ...) 
{
    get.g <- function(x, data, ...) {
        g.stats.glob(data.frame(x, data), ...)$g.stats
    }
#    nobs <- length(level)
    perm.stat <- vector(length = nperm)
    perm.stat[nperm] <- get.g(level, data, ...)
    for (i in 1:(nperm - 1)) {
        perm.stat[i] <- get.g(sample(level), data[, ], ...)
    }
    list(g.star = perm.stat, p.val = sum(perm.stat >= perm.stat[nperm])/nperm)
}
"test.within" <-
function (data = data, within, test.lev, nperm = 100, ...) 
{
    get.g <- function(x, data, ...) {
        g.stats.glob(data.frame(x, data), ...)$g.stats
    }
    within<-as.integer(factor(within))
    test.lev<-as.integer(factor(paste(within,"X",test.lev,sep="")))
    x<-order(within,test.lev)
    data<-data.frame(data[x,])
    within<-within[x]
    test.lev<-test.lev[x]
#    tlev<-paste(within,test.lev,sep="")
#    tlev<-rep(1:dim(table(tlev)),c(table(tlev)))
#    nobs <- length(test.lev)
    perm.stat <- vector(length = nperm)
    perm.stat[nperm] <- get.g(test.lev, data, ...)
    for (i in 1:(nperm - 1)) {
        perm.stat[i] <- get.g(test.lev, data[samp.within(within), 
            ], ...)
    }
    list(g.star = perm.stat, p.val = sum(perm.stat >= perm.stat[nperm])/nperm)
}
