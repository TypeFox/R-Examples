###############################################################
################ Fts class definitions ####################
###############################################################
###############################################################

## layout of an fts object is
## a matrix with a dates attribute attached
fts <- function(index,data) {
    if(missing(index)) stop("missing index")
    if(missing(data)) stop("missing data")
    stopifnot(length(index)==NROW(data))
    structure(as.matrix(data),index=index,class=c("fts","zoo"))
}

as.fts <- function(x) {
    UseMethod("as.fts")
}

as.fts.default <- function(x) {
    if(!is.null(dim(x))) {
        dts <- rownames(x)
    } else {
        dts <- names(x)
    }
    fts(index=dts,data=x)
}

as.fts.data.frame <- function(x) {
    cnames <-  colnames(x)
    ans <- fts(index=x[,"asofdate"],
               data=as.matrix(x[,-match("asofdate",cnames)]))
    colnames(ans) <- cnames[-1]
    ans
}

as.fts.zoo <- function(x) fts(index=attr(x,"index"),data=as.matrix(coredata(x)))

as.matrix.fts <- function(x, ...)
    structure(matrix(as.numeric(x),nrow=nrow(x),ncol=ncol(x)),dimnames=list(format(index(x),"%Y%m%d"),colnames(x)))

as.dataframe.fts <- function(x, ...) {
    ans <- data.frame(x)
    rownames(ans) <- format(index(x),"%Y%m%d")
    ans
}

## create an fts object given dates and column names
template.fts <- function(index,cnames) {
    ans <- fts(index=index,
               data=matrix(nrow=length(index),ncol=length(cnames)))
    colnames(ans) <- cnames
    ans
}

Ops.fts <- function (e1, e2) {

    if(missing(e2)) {
        .Class <- "matrix"
        NextMethod()
    } else {
        c.e1 <- class(e1)
        c.e2 <- class(e2)

        if("fts" %in% c.e1 && "fts" %in% c.e2) {
            nce1 <- NCOL(e1)
            nce2 <- NCOL(e2)

            if(nce1!=nce2 && nce1!=1 && nce2!=1) {
                stop("Ops.fts: non conformable data.")
            }
            stopifnot(all.equal(class(index(e1)),class(index(e2))))
            i.dates <- intersect(index(e1),index(e2))
            class(i.dates) <- class(index(e1))

            ## if there is an intersection, the do the Op
            if(length(i.dates)) {
                e1 <- e1[i.dates,]
                e2 <- e2[i.dates,]

                if(nce1==1 && nce2!=1) {
                    e1 <- rep(e1,nce2)
                } else if(nce1!=1 && nce2==1) {
                    e2 <- rep(e2,nce1)
                }

                .Class <- "matrix"
                ans <- NextMethod()
                attr(ans,"index") <- i.dates
                class(ans) <- c("fts","zoo")
            } else {
                ## no matching dates, return NULL
                ans <- NULL
            }
        } else {
            .Class <- "matrix"
            ans <- NextMethod()
            if("fts" %in% c.e1) {
                ans.dates <- attr(e1,"index")
            } else {
                ans.dates <- attr(e2,"index")
            }
            attr(ans,"index") <- ans.dates
            class(ans) <- c("fts","zoo")
        }
        ans
    }
}

print.fts <- function(x, ...) {
    cnms <- colnames(x)
    dnms <- list(format(index(x)),cnms)
    dims <- dim(x)
    attributes(x) <- NULL
    dim(x) <- dims
    dimnames(x) <- dnms
    print(x,...)
    invisible(x)
}

"[.fts" <- function(x,i,j,...,drop=FALSE) {

    if(missing(i)) i <- 1:nrow(x)
    if(missing(j)) j <- 1:ncol(x)

    ##
    ##if(typeof(i)=="character") {
    ##    i <- guess.index(i,index(x))
    ##}

    ## if we have dates, then use them
    if(any(c("Date","POSIXct") %in% class(i))) {
        i <- match(i,index(x))
    }

    ## FIXME: not sure about this hack...
    ##i <- i[!is.na(i)]
    ##j <- j[!is.na(j)]

    ans.dates <- index(x)[i]
    ans <- unclass(x)[i,j,...,drop=drop]
    attr(ans,"index") <- ans.dates
    class(ans) <- c("fts","zoo")

    ans
}

"[<-.fts" <- function(x, i = TRUE, j = TRUE, ..., value) {
    ## if we have dates, then use them
    if(!missing(i)) {
        if(any(c("Date","POSIXct") %in% class(i))) {
            i <- match(i,index(x))
        } else if (mode(i)=="logical") {
            i <- which(i)
        }
    }
    .Class <- "matrix"
    NextMethod()
}

## apply a function to an fts row
row.apply <- function(x,FUN,...) {
    ans <- apply(x,1,FUN,...)

    if(!is.null(dim(ans))) {
        ans <- t(ans)
        rownames(ans) <- NULL
    }

    fts(index=index(x),data=ans)        
}

row.any <- function(x) apply(x,1,any)
row.all <- function(x) apply(x,1,all)
col.any <- function(x) apply(x,2,any)
col.all <- function(x) apply(x,2,all)

## apply a function to an fts column
column.apply <- function(x,FUN,...) {
    apply(x,2,FUN,...)
}

remove.na.rows <- function(x) x[!row.any(is.na(x)),]

remove.all.na.rows <- function(x) x[!row.all(is.na(x)),]

as.data.frame.fts <- function(x,row.names = NULL, optional = FALSE, check.names = TRUE,...) {
    data.frame(asofdate=index(x),as.data.frame.matrix(unclass(x),optional=optional),check.names=check.names)
}

as.matrix.fts <- function(x, ...) {
    rownames(x) <- format(index(x))
    attr(x,"index") <- NULL
    class(x) <- "matrix"
    x
}

rbind.fts <- function(...) {
    x <- list(...)

    ## check for same number of cols
    if(length(unique(unlist(lapply(x,ncol))))!=1) {
        stop("unequal number of cols in arguments.")
    }

    ## unclass data
    ans.unclass <- lapply(x,unclass)

    ## do rbind
    ans.data <- do.call(rbind,ans.unclass)

    ## get dates
    ans.dates <- do.call(c,lapply(x,index))

    ## sort it
    new.order <- order(ans.dates)
    ## must do data before we sort dates
    ans.data <- ans.data[new.order,,drop=F]
    ans.dates <- ans.dates[new.order]
    fts(index=ans.dates,data=ans.data)
}

cbind.fts <- function(...) {
    x <- list(...)

    ans.index <- sort(unique(unlist(lapply(x,index))))
    class(ans.index) <- class(index(x[[1]]))
    nrow.ans <- length(ans.index)
    ncol.list <- lapply(x,ncol)
    ncol.ans <- sum(unlist(ncol.list))

    ans <- matrix(NA,nrow=nrow.ans,ncol=ncol.ans)

    ## starting col
    start.col <- 1

    ## fill the answer
    for(i in 1:length(x)) {
        ## guard against null contracts from LIM
        if(!is.null(x[[i]])) {
            end.col <- start.col + ncol.list[[i]] - 1
            ans[match(index(x[[i]]),ans.index),start.col:end.col] <- x[[i]]
            start.col <- end.col + 1
        }
    }

    ## set attributes
    ## fix blank cnames
    cnames.list <- lapply(x,colnames)
    colnames(ans) <- fix.cnames(cnames.list)
    attr(ans,"index") <- ans.index
    class(ans) <- c("fts","zoo")
    ans
}

fix.cnames <- function(cnames.list) {
    null.cnames <- unlist(lapply(cnames.list,is.null))
    for(i in 1:length(cnames.list)) {
        if(null.cnames[i]) {
            cnames.list[[i]] <- ""
        }
    }
    unlist(cnames.list)
}

trim <- function(x,trim.dates) {
    new.dates <- sort(intersect(index(x),trim.dates))
    class(new.dates) <- class(index(x))
    x[new.dates,]
}

write.csv.fts <- function(x, file, ...) {
    write.csv(as.data.frame(x), file,row.names=FALSE, ...)
}

read.csv.fts <- function(file, date.column=1, date.format="%Y-%m-%d",date.convert.fun=as.Date,...) {
    fts.data <- read.csv(file,...)

    if(mode(date.column)=="character")
        date.column <- match(date.column,colnames(fts.data))

    fts(index=date.convert.fun(strptime(fts.data[,date.column],date.format)),
        data=as.matrix(fts.data[, -date.column, drop=F]))
}

cumsum.fts <- function(x) fts(index=index(x),data=apply(x,2,cumsum))
cumprod.fts <- function(x) fts(index=index(x),data=apply(x,2,cumprod))
cummax.fts <- function(x) fts(index=index(x),data=apply(x,2,cummax))
cummin.fts <- function(x) fts(index=index(x),data=apply(x,2,cummin))

###############################################################
################### Date Functions ############################
###############################################################
###############################################################
## return the dates just when the fts is true
event.dates <- function(x) {
    stopifnot(ncol(x)==1)
    index(x)[as.logical(x)&!is.na(x)]
}

## find date intersection of all tseries
intersect.all <- function(...) {
    x <- list(...)
    dts <- lapply(x,index)
    ans <- dts[[1]]

    if(length(dts) > 1) {
        for(i in 2:length(x)) {
            ans <- intersect(ans,dts[[i]])
        }
    }

    class(ans) <- class(index(x[[1]]))
    ans
}


###############################################################
################ Calls to External Library ####################
###############################################################
###############################################################


moving.mean <- function(x,periods) {
    .Call("movingMean",x,as.integer(periods),PACKAGE="fts")
}

moving.sum <- function(x,periods) {
    .Call("movingSum",x,as.integer(periods),PACKAGE="fts")
}

moving.product <- function(x,periods) {
    .Call("movingProduct",x,as.integer(periods),PACKAGE="fts")
}

moving.max <- function(x,periods) {
    .Call("movingMax",x,as.integer(periods),PACKAGE="fts")
}

moving.min <- function(x,periods) {
    .Call("movingMin",x,as.integer(periods),PACKAGE="fts")
}

moving.rank <- function(x,periods) {
    .Call("movingRank",x,as.integer(periods),PACKAGE="fts")
}

moving.sd <- function(x,periods) {
    .Call("movingStdev",x,as.integer(periods),PACKAGE="fts")
}

expanding.max <- function(x) {
    .Call("expandingMax",x,PACKAGE="fts")
}

expanding.min <- function(x) {
    .Call("expandingMin",x,PACKAGE="fts")
}


moving.cor <- function(x,y,periods) {
    .Call("movingCor", x, y, as.integer(periods),PACKAGE="fts")
}

moving.cov <- function(x,y,periods) {
    .Call("movingCov", x, y, as.integer(periods),PACKAGE="fts")
}


since.na <- function(x) {
    .Call("sinceNA",x,PACKAGE="fts")
}

lag.fts <- function(x, k, ...) {
    stopifnot(k > -1)
    ans <- .Call("lag", x, as.integer(k),PACKAGE="fts")
    ans
}

lead <- function(x, k, ...) {
    UseMethod("lead")
}

lead.fts <- function(x, k, ...) {
    if(k < 0) stop("only positive values of k are allowed")
    .Call("lead",x ,as.integer(k),PACKAGE="fts")
}

diff.fts <- function(x, k, ...) {
    stopifnot(k > 0)
    ans <- .Call("diff", x, as.integer(k),PACKAGE="fts")
    attr(ans,"ticker") <- attr(x,"ticker")
    ans
}

fill.fwd <- function(x) {
    .Call("fillForward",x,PACKAGE="fts")
}

fill.bwd <- function(x) {
    .Call("fillBackward",x,PACKAGE="fts")
}

fill.value <- function(x,value) {
    .Call("fillValue",x,value,PACKAGE="fts")
}

pad <- function(x,pad.dates) {
    stopifnot(all.equal(class(index(x)),class(pad.dates)))
    if(storage.mode(pad.dates)!=storage.mode(index(x))) storage.mode(pad.dates) <- storage.mode(index(x))

    .Call("pad",x,pad.dates,PACKAGE="fts")
}

monthly.sum <- function(x) {
    .Call("monthlySum",x,PACKAGE="fts")
}

to.yearly <- function(x) {
    .Call("toYearly",x,PACKAGE="fts")
}

to.quarterly <- function(x) {
    .Call("toQuarterly",x,PACKAGE="fts")
}

to.monthly <- function(x) {
    .Call("toMonthly",x,PACKAGE="fts")
}

to.weekly <- function(x) {
    .Call("toWeekly",x,PACKAGE="fts")
}

to.daily <- function(x) {
    .Call("toDaily",x,PACKAGE="fts")
}

to.hourly <- function(x) {
    .Call("toHourly",x,PACKAGE="fts")
}

to.minute <- function(x) {
    .Call("toMinute",x,PACKAGE="fts")
}

to.second <- function(x) {
    .Call("toSecond",x,PACKAGE="fts")
}

to.day.of.week <- function(x,day.of.week,beginning.of.period=TRUE) {
    dts <- index(x)
    stopifnot(class(dts)=="Date")
    end.date <- dts[nrow(x)]
    end.date <- end.date + 6
    pad.days <- seq(from=dts[1],to=end.date,by="days")
    ## seq changes storage mode
    storage.mode(pad.days) <- storage.mode(dts)
    pad.days <- pad.days[as.POSIXlt(pad.days)$wday==day.of.week]
    x.filled <- fill.fwd(pad(x,pad.days))
    ans <- x.filled[as.POSIXlt(index(x.filled))$wday == day.of.week,]
    if(beginning.of.period) {
        index(ans) <- index(ans) - 6
    }
    ans
}

###############################################################
############ Plotting functions for Fts Objects ###########
###############################################################
###############################################################

plot.fts <- function(x, type="l", xlab="Date", ylab=substitute(x), ...) {
    x.range <- c(min(index(x)), max(index(x)))
    y.range <- c(min(x, na.rm=TRUE), max(x, na.rm=TRUE))
    plot(x.range, y.range, type="n", xlab=xlab, ylab=ylab, ...)

    params <- lapply(list(type = type, ...), rep, length.out = ncol(x))
    for (i in seq(ncol(x))) {
        args <- c(list(xy.coords(index(x), x[,i])),
                  lapply(params, `[`, i))
        do.call(plot.xy, args)
    }
}

## drop out rows that do not have the required number of observations
filter.min.obs <- function(x,obs.required) {
    obs <- apply(!is.na(x),1,sum)
    x[obs >= obs.required,]
}

###############################################################
################### Technical Analysis  #######################
###############################################################
###############################################################
ema <- function(x,periods) {
    .Call("ema",x,as.integer(periods),PACKAGE="fts")
}

rsi <- function(x,periods) {
    x.up <- x
    x.up[x<0] <- 0
    x.down <- x
    x.down[x>0] <- 0

    x.avg.up <- ema(x.up, periods)
    x.avg.down <- ema(x.down, periods)

    100 - 100/(1 - x.avg.up/x.avg.down)
}

year <- function(x) {
    as.POSIXlt(index(x))$year+1900
}

month <- function(x) {
    as.POSIXlt(index(x))$mon+1
}

mday <- function(x) {
    as.POSIXlt(index(x))$mday
}

wday <- function(x) {
    as.POSIXlt(index(x))$wday
}

ma.crossover.down <- function(x,n) {
    has.close <- !is.null(colnames(x)) && "close" %in% colnames(x)
    stopifnot(has.close || ncol(x) == 1)
    col <- ifelse(has.close,"close",1)

    x <- x[,col]

    xma <- moving.mean(x,n)

    x < xma & lag(x,1) > lag(xma,1)
}

ma.crossover.up <- function(x,n) {
    has.close <- !is.null(colnames(x)) && "close" %in% colnames(x)
    stopifnot(has.close || ncol(x) == 1)
    col <- ifelse(has.close,"close",1)

    x <- x[,col]

    xma <- moving.mean(x,n)

    x > xma & lag(x,1) < lag(xma,1)
}

ma.crossover <- function(x,n) {
    ma.crossover.up(x,n) - ma.crossover.down(x,n)
}

lower.low <- function(x) {
    stopifnot("low" %in% colnames(x))
    xl <- x[,"low"]
    xl < lag(xl,1)
}

higher.high <- function(x) {
    stopifnot("high" %in% colnames(x))
    xh <- x[,"high"]
    xh > lag(xh,1)
}

repeated <- function(x, times) {
    stopifnot(ncol(x)==1 && mode(x)=="logical")
    moving.sum(x,times)==as.integer(times)
}

new.low <- function(x,n) {
    stopifnot("low" %in% colnames(x))
    moving.rank(x[,"low"],n)==1
}

new.high <- function(x,n) {
    stopifnot("high" %in% colnames(x))
    moving.rank(x[,"high"],n)==as.integer(n)
}

above.ma <- function(x,n) {
    stopifnot("close" %in% colnames(x))
    xc <- x[,"close"]
    xc > moving.mean(xc,n)
}

below.ma <- function(x,n) {
    stopifnot("close" %in% colnames(x))
    xc <- x[,"close"]
    xc < moving.mean(xc,n)
}

ma.d <- function(x,n) {
    above.ma(x,n) - below.ma(x,n)
}

higher.low <- function(x) {
    stopifnot("low" %in% colnames(x))
    xl <- x[,"low"]
    xl > lag(xl,1)
}

lower.high <- function(x) {
    stopifnot("high" %in% colnames(x))
    xh <- x[,"high"]
    xh < lag(xh,1)
}

up <- function(x) {
    stopifnot("close" %in% colnames(x))
    xc <- x[,"close"]
    xc > lag(xc,1)
}

down <- function(x) {
    stopifnot("close" %in% colnames(x))
    xc <- x[,"close"]
    xc < lag(xc,1)
}

pct.chg <- function(x) {
    stopifnot("close" %in% colnames(x))
    xc <- x[,"close"]
    diff(xc,1) / lag(xc,1) * 100
}

inside.day <- function(x) {
    stopifnot(all(c("high","low") %in% colnames(x)))

    xh <- x[,"high"]
    xl <- x[,"low"]

    xh < lag(xh,1) & xl > lag(xl,1)
}

inside.day.up <- function(x) {
    inside.day(x) & up(x)
}

inside.day.down <- function(x) {
    inside.day(x) & down(x)
}

inside.day.direction <- function(x) {
    inside.day.up(x) - inside.day.down(x)
}

outside.day <- function(x) {
    stopifnot(all(c("high","low") %in% colnames(x)))

    xh <- x[,"high"]
    xl <- x[,"low"]

    xh > lag(xh,1) & xl < lag(xl,1)
}

outside.day.up <- function(x) {
    outside.day(x) & up(x)
}

outside.day.down <- function(x) {
    outside.day(x) & down(x)
}

outside.day.direction <- function(x) {
    outside.day.up(x) - outside.day.down(x)
}

hl.oc.ratio <- function(x) {
    stopifnot(all(c("open","high","low","close") %in% colnames(x)))
    abs(x[,"close"] - x[,"open"]) / (x[,"high"] - x[,"low"])
}

gap.up <- function(x) {
    stopifnot(all(c("open","high","low","close") %in% colnames(x)))
    xh <- x[,"high"]
    xo <- x[,"open"]
    xo > lag(xh,1)
}

gap.down <- function(x) {
    stopifnot(all(c("open","high","low","close") %in% colnames(x)))
    xl <- x[,"low"]
    xo <- x[,"open"]
    xo < lag(xl,1)
}

gap.up.continue <- function(x) {
    gap.up(x) & up(x)
}

gap.down.continue <- function(x) {
    gap.down(x) & down(x)
}

gap.continue <- function(x) {
    gap.up.continue(x) - gap.down.continue(x)
}

gap.up.reverse <- function(x) {
    gap.up(x) & down(x)
}

gap.down.reverse <- function(x) {
    gap.down(x) & up(x)
}

gap.reverse <- function(x) {
    gap.up.reverse(x) - gap.down.reverse(x)
}

gap.direction <- function(x) {
    gap.up(x) - gap.down(x)
}

rsi.crossover.up <- function(x,periods,thresh) {
    stopifnot("close" %in% colnames(x))
    xc <- x[,"close"]
    xd <- diff(xc,1)
    x.rsi <- rsi(xd,periods)
    x.rsi <- x.rsi[!is.na(x.rsi),]
    x.rsi < thresh & lag(x.rsi,1) > thresh
}

rsi.crossover.down <- function(x, periods, thresh) {
    stopifnot("close" %in% colnames(x))
    xc <- x[,"close"]
    xd <- diff(xc,1)
    x.rsi <- rsi(xd,periods)
    x.rsi <- x.rsi[!is.na(x.rsi),]
    x.rsi > thresh & lag(x.rsi,1) < thresh
}

rsi.crossover <- function(x,periods, thresh) {
    rsi.crossover.up(x,periods,thresh) - rsi.crossover.down(x,periods,thresh)
}

ma.distance <- function(x, periods) {
    stopifnot("close" %in% colnames(x))
    xc <- x[,"close"]
    xma <- moving.mean(xc,periods)
    (xc - xma) / xc * 100
}

trend.day.up <- function(x,thresh=0.2) {
    stopifnot(all(c("open","high","low","close") %in% colnames(x)))
    xo <- x[,"open"]
    xh <- x[,"high"]
    xl <- x[,"low"]
    xc <- x[,"close"]
    day.range <- xh - xl
    open.low <- (xo - xl) / day.range
    high.close <- (xh - xc) / day.range
    open.low <= thresh & high.close <= thresh
}

trend.day.down <- function(x,thresh=0.2) {
    stopifnot(all(c("open","high","low","close") %in% colnames(x)))
    xo <- x[,"open"]
    xh <- x[,"high"]
    xl <- x[,"low"]
    xc <- x[,"close"]
    day.range <- xh - xl
    open.high <- (xh - xo) / day.range
    low.close <- (xc - xl) / day.range
    open.high <= thresh & low.close <= thresh
}

trend.day <- function(x,thresh=.2) {
    trend.day.up(x,thresh) - trend.day.down(x,thresh)
}

cor.by.row <- function(x,y) {
    i.dts <- sort(intersect(index(x),index(y)))
    class(i.dts) <- class(index(x))
    ans <- template.fts(i.dts,"cor")

    for(i in 1:length(i.dts)) {
        ans[i,] <- cor(as.vector(x[i.dts[i],]),as.vector(y[i.dts[i],]))
    }
    ans
}
