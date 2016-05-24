##    itspkg.r: a time series package in R
##    Copyright (C) 2003 to 2004 Commerzbank Securities
##    Copyright (C) 2004 to present Whit Armstrong


itsState <- new.env()
assign(x=".itsformat", value="%Y-%m-%d" , env=itsState)
setClass("its",representation("matrix",dates="POSIXt"))

##-Methods-
##names-method-------------------------------------------------------
namesIts <- function(x) {
    return(dimnames(x@.Data)[[2]])
}
setMethod("names",signature(x="its"),namesIts)
##"names<-"-method---------------------------------------------------
"names<-Its" <- function(x,value) {
    result <- x
    dimnames(result)[[2]] <- value
    return(result)
}

setMethod("names<-",signature(x="its",value="character"),get("names<-Its"))
##dates-method-------------------------------------------------------
if(!isGeneric("dates")) {setGeneric("dates", function(x,...) standardGeneric("dates"))}
datesIts <- function(x) {
    return(x@dates)
}
setMethod("dates",signature(x="its"),datesIts)
##"dates<-"-method---------------------------------------------------
if(!isGeneric("dates<-")) {setGeneric("dates<-", function(x,value) standardGeneric("dates<-"))}
"dates<-Its" <- function(x,value) {
    if(!is(value,"POSIXt")) stop("dates should be in POSIX format")
    value <- as.POSIXct(value)
    result <- its(x@.Data,value)
    return(result)
}
setMethod("dates<-",signature(x="its",value="POSIXt"),get("dates<-Its"))
##core-method-------------------------------------------------------
if(!isGeneric("core")) {setGeneric("core", function(x) standardGeneric("core"))}
coreIts <- function(x) {
    return(x@.Data)
}
setMethod("core",signature(x="its"),coreIts)
##"core<-"-method---------------------------------------------------
if(!isGeneric("core<-")) {setGeneric("core<-", function(x,value) standardGeneric("core<-"))}
"core<-Its" <- function(x,value) {
    result <- its(value,dates(x))
    names(result) <- names(x)
    return(result)
}
setMethod("core<-",signature(x="its",value="matrix"),get("core<-Its"))
##arith-methods------------------------------------------------------
its.its.opp <- function(e1,e2) {
    ## take intersection of dates
    i.dates <- sort(intersect(dates(e1),dates(e2)))
    class(i.dates) <- c("POSIXt","POSIXct")

    ## add the data, taking the subset of the core for which the dates match
    ans <- callGeneric(e1[dates=i.dates,]@.Data,e2[dates=i.dates,]@.Data)

    ## make a new its w/ the ans and the dates intersection
    return(its(ans,i.dates))
}
setMethod("Arith",signature(e1="its",e2="its"),its.its.opp)

its.numeric.opp <- function(e1,e2) {
    ans <- callGeneric(e1@.Data,e2)
    return(its(ans,dates(e1)))
}
setMethod("Arith",signature("its", "numeric"),its.numeric.opp)

numeric.its.opp <- function(e1,e2) {
    ans <- callGeneric(e1,e2@.Data)
    return(its(ans,dates(e2)))
}
setMethod("Arith",signature("numeric","its"),numeric.its.opp)
##plot-method--------------------------------------------------------
if(!isGeneric("plot")) setGeneric("plot", useAsDefault=plot)

plotIts <- function(x,y,colvec=1:ncol(x),type="l",ltypvec=1,lwdvec=1,
                    leg=FALSE,yrange,format,at,interp=c("linear","none"),lab=FALSE,...)
{
    if(missing(yrange)){ylim <- range(x,na.rm=TRUE)} else {ylim <- yrange}
    interp <- match.arg(interp)
    firstp <- TRUE
    xdates <- x@dates
    n <- dim(x)[1]
    m <- dim(x)[2]
    ##make line control parameters correct length
    colveclong <- rep(colvec,length.out=m)
    ltypveclong <- rep(ltypvec,length.out=m)
    lwdveclong <- rep(lwdvec,length.out=m)
    for(i in 1:m)
    {
        if(interp=="linear")
        {
            vpoints <- c(1,which(!is.na(x[,i])),n)
            xxx <- x[,i]
        } else
    {
        vpoints <- 1:n
        xxx <- expandIts(x[,i])
    }
        for (j in 1:ncol(xxx))
        {
            if(!firstp){par(new=TRUE)}else {firstp <- FALSE}
            plot(x=xdates[vpoints],
                 y=xxx[vpoints,j],
                 type=type,
                 col=colveclong[i],
                 ylim=ylim,
                 lty=ltypveclong[i],
                 lwd=lwdveclong[i],
                 xaxt="n",
                 ...)
        }
    }
    if(lab)
    {
        labcurve(curves=gencurves(x),
                 labels=dimnames(x)[[2]],
                 col=colveclong,
                 cex=.8)
    } else if(leg)
    {
        labcurve(curves=gencurves(x),
                 labels=dimnames(x)[[2]],
                 col=colveclong,
                 lty=ltypveclong[i],
                 lwd=lwdveclong[i],
                 keys=rep(" ",ncol(x)), ##letters[1:ncol(x)],
                 cex=.8)
    }
    grid()
    axis.POSIXct(x=xdates[vpoints],side=1,at=at,format=format)
}
setMethod("plot",signature(x="its",y="missing"),plotIts)

##print-method-------------------------------------------------------
printIts <- function(x,...){
    print(x@.Data,...)
}
setMethod("print",signature(x="its"),printIts)

##start-method-------------------------------------------------------
startIts <- function(x,format=its.format(),...) {
    return(format(x@dates[1],format=format,...))
}
setMethod("start",signature(x="its"),startIts)

##end-method--------------------------------------------------------
endIts <- function(x,format=its.format(),...) {
    return(format(x@dates[length(x@dates)],format=format,...))
}
setMethod("end",signature(x="its"),endIts)

##summary-method-----------------------------------------------------
summaryIts <- function(object,...) {
    r1 <- apply(object,2,min,na.rm=TRUE)
    r2a <- apply(object,2,quantile,probs=.25,na.rm=TRUE)
    r2b <- apply(object,2,quantile,probs=.5,na.rm=TRUE)
    r3 <- apply(object,2,mean,na.rm=TRUE)
    r4 <- apply(object,2,quantile,probs=.75,na.rm=TRUE)
    r5 <- apply(object,2,max,na.rm=TRUE)
    r6 <- apply(is.na(object),2,sum,na.rm=TRUE)
    r7 <- rep(nrow(object),ncol(object))-r6
    r8 <- sqrt(apply(object,2,var,na.rm=TRUE))

    mysum <- rbind(r1,r2a,r2b,r3,r4,r5,r6,r7,r8)
    dimnames(mysum)[[1]] <- c("Min.","1st Qu.","Median","Mean","3rd Qu.","Max.","NA's","non-NA's","s.d.")
    dimnames(mysum)[[2]] <- names(object)
    mysum
}

setMethod("summary",signature(object="its"),summaryIts)

##cumsum-method------------------------------------------------------
cumsumIts <- function(x) {
    y <- x
    for(j in 1:ncol(x)) {y[,j] <- cumsum(as.numeric(x[,j]))}
    return(y)
}
setMethod("cumsum",signature(x="its"),cumsumIts)

##diff-method--------------------------------------------------------
diffIts <- function(x,lag=1) {
    i <- NULL
    if((lag+1)<=length(x@dates)) i <- (lag+1):length(x@dates)
    y <- its(diff(x@.Data,lag=lag,drop=FALSE),dates=x@dates[i])
    return(y)
}
setMethod("diff",signature(x="its"),diffIts)

##union-method-------------------------------------------------------
unionIts <- function(x,y)
{
    if(!is.null(x)&!is.null(y))
    {
        dates1                  <- x@dates
        dates2                  <- y@dates
        alldates                <- sort(union(dates1,dates2))
        allnames                <- c(dimnames(x)[[2]],dimnames(y)[[2]])
        class(alldates)         <- class(x@dates)
        m1                      <- ncol(x)
        m2                      <- ncol(y)
        n                       <- length(alldates)
        m                       <- m1+m2
        united                  <- matrix(NA,nrow=n,ncol=m)
        if(m1>0) united[match(dates1,alldates),1:m1] <- x
        if(m2>0) united[match(dates2,alldates),(m1+1):m] <- y
        result <- its(united,dates=alldates,names=allnames)
    }
    if(is.null(x)) {result <- y}
    if(is.null(y)) {result <- x}
    return(result)
}
setMethod("union",signature(x="its",y="its"),unionIts)
setMethod("union",signature(x="its",y="NULL"),unionIts)
setMethod("union",signature(x="NULL",y="its"),unionIts)

##intersect-method---------------------------------------------------
intersectIts <- function(x,y) {
  if( !is.null(x) & !is.null(y) ) {
    dates1                  <- x@dates
    dates2                  <- y@dates
    alldates                <- sort(intersect(dates1,dates2))
    class(alldates)         <- class(x@dates)
    allnames                <- c(dimnames(x)[[2]],dimnames(y)[[2]])
    m1                      <- dim(x)[2]
    m2                      <- dim(y)[2]
    n                       <- length(alldates)
    m                       <- m1+m2
    united                  <- matrix(NA,n,m)
    drows1                  <- sort(match(dates1,alldates))
    drows2                  <- sort(match(dates2,alldates))
    srows1                  <- sort(match(alldates,dates1))
    srows2                  <- sort(match(alldates,dates2))
    united[drows1,1:m1]     <- x[srows1,,drop=FALSE]
    united[drows2,(m1+1):m] <- y[srows2,,drop=FALSE]
    result <- its(united,dates=alldates,names=allnames)
  }
  if(is.null(x)) {result <- y}
  if(is.null(y)) {result <- x}
  return(result)
}
setMethod("intersect",signature(x="its",y="its"),intersectIts)
setMethod("intersect",signature(x="its",y="NULL"),intersectIts)
setMethod("intersect",signature(x="NULL",y="its"),intersectIts)

##as-method--------------------------------------------------------
as.data.frameIts <- function(from) {
    data.frame(core(from))
}
setAs(from="its",to="data.frame",as.data.frameIts)

##-Functions-
##readcsvIts-function------------------------------------------------
readcsvIts <- function(filename,informat=its.format(),outformat=its.format(),tz="",usetz=FALSE,header=TRUE,...)
{
    mydata <- read.csv(filename,header=header,...)
    n <- dim(mydata)[1]
    m <- dim(mydata)[2]
    datamat <- as.numeric(as.matrix((mydata)[,2:m,drop=FALSE]))
    dim(datamat) <- c(n,(m-1))
    dimnames(datamat) <- list(dimnames(mydata)[[1]],dimnames(mydata)[[2]][2:m])
    dimnames(datamat)[[1]]  <- format(strptime(as.character(as.vector(mydata[,1])),informat),
                                      format=outformat,tz=tz,usetz=usetz)
    return(datamat)
}
##writecsvIts-function-----------------------------------------------
writecsvIts <- function(x,filename,format=its.format(),tz="",usetz=FALSE,col.names=NA,sep=",",split=FALSE,...)
{
    if (!inherits(x, "its")) stop("function is only valid for objects of class 'its'")
    dimnames(x)[[1]] <- format(x@dates,format=format,tz=tz,usetz=usetz)
    y <- data.frame(x@.Data)
    dimnames(y) <- dimnames(x)
    if(split & ncol(x)>255)
    {
        jstart <- 1
        jend <- 255
        j <- 0
        while(jstart<ncol(x))
        {
            j <- j+1
            fnam <- paste(strsplit(filename,"\\.")[[1]],collapse=paste(j,".",sep=""))
            write.table(y[,jstart:jend],file=fnam,col.names=col.names,sep=sep,...)
            jstart <- jend+1
            jend <- min(c(jstart+254,ncol(x)))
        }
    } else
{
    mydata <- write.table(y,file=filename,col.names=col.names,sep=sep,...)
}
}

##accrueIts-function-------------------------------------------------
accrueIts <- function(x,daysperyear=365)
{
    if (!inherits(x, "its")) stop("function is only valid for objects of class 'its'")
    n <- nrow(x)
    accrued.val <- x[-n,]*diff(as.numeric(x@dates))/(daysperyear*24*60*60)
    accrued.dates <- x@dates[-1]
    return(its(accrued.val,dates=accrued.dates))
}
##its-function-------------------------------------------------------
its <- function(x,
                dates=as.POSIXct(x=strptime(dimnames(x)[[1]],format=its.format())),
                names=dimnames(x)[[2]],format=its.format(),...)
{

    if(!is(dates,"POSIXt")) stop("dates should be in POSIX format")

    dates <- as.POSIXct(dates)

    ## fix identical bug
    if(is.null(attr(dates,"tzone"))) attr(dates,"tzone") <- ""

    if(is.null(dim(x))){dim(x) <- c(length(x),1)}
    x <- addDimnames(x)
    if(!(nrow(x)==length(dates))) {stop("dates length must match matrix nrows")}
    if(!(ncol(x)==length(names))) {stop("names length must match matrix ncols")}
    dimnames(x)[[1]] <- format(dates,format=format,...)
    dimnames(x)[[2]] <- names
    return(new("its",x,dates=dates))
}
##is.its-function----------------------------------------------------
is.its <- function(object)
{
    return(inherits(object,"its") && validIts(object))
}


##as.its-method-------------------------------------------------------
as.its <- function(x,...) { UseMethod("as.its") }

## as.its-fucntion--------------------------------------------
as.its.default <- function(x,...)
{
    dates <- as.vector(x[,1])
    class(dates) <- c("POSIXt","POSIXct")
    return(its(x=x[,-1],dates=dates,...))
}

## as.its.zoo-fucntion--------------------------------------------
## for converting an its object into a zoo object
## contributed by the zoo team
as.its.zoo <- function(x,...)
{
    stopifnot(require(its))
    index <- attr(x, "index")
    stopifnot(inherits(index, "POSIXct"))
    attr(x, "index") <- NULL
    its(unclass(x), index)
}

##lagIts-function----------------------------------------------------
lagIts <- function(x,k=1)
{

    if (!inherits(x, "its")) stop("function is only valid for objects of class 'its'")

    lagmat <- core(x)*NA

    dimnames(lagmat)[[2]] <- paste(dimnames(lagmat)[[2]],"lag",k)

    n <- dim(x)[1]

    if(k>0) {
        lagmat[(k+1):n,] <- x[1:(n-k),]
    }  else {
        lagmat[1:(n+k),] <- x[(1-k):n,]
    }

    y <- its(lagmat,dates=x@dates)

    return(y)
}

##lagdistIts-function------------------------------------------------
lagdistIts <- function(x,kmin,kmax)
{
    if (!inherits(x, "its")) stop("function is only valid for objects of class 'its'")
    result <- lagIts(x,kmin)
    if(kmax>kmin)
    {
        for(j in (kmin+1):kmax)
        {
            result <- intersectIts(result,lagIts(x,j))
        }
    }
    return(result)
}

##alignedIts-function---------------------------------------------
alignedIts <- function(obj1,obj2,print=FALSE) {
    if( !inherits(obj1, "its") & inherits(obj2, "its") )
      stop("function is only valid for objects of class 'its'")

    ##takes the intersection of the dates and extracts same dates for both
    mat <- intersectIts(obj1,obj2)
    obj1a <- mat[,1:ncol(obj1),drop=FALSE]
    obj2a <- mat[,(ncol(obj1)+1):ncol(mat),drop=FALSE]
    if(print) {
      print(paste("inputs number of rows",nrow(obj1),nrow(obj2),"; output number of rows",nrow(mat)))
      print(paste("inputs number of cols",ncol(obj1),ncol(obj2),"; output number of cols",ncol(obj1a),ncol(obj2a)))
    }
    list(obj1a,obj2a)
}

##appendIts-function-------------------------------------------------
appendIts <- function(obj1,obj2,but=TRUE,matchnames=TRUE)
{
    if(is.null(obj1)&inherits(obj2, "its")) return(obj2)
    if(is.null(obj2)&inherits(obj1, "its")) return(obj1)
    if (!inherits(obj1, "its")&inherits(obj2, "its")) stop("function is only valid for objects of class 'its'")

    overlap <- overlapsIts(obj1,obj2)

    if(overlap & but) stop("overlap not allowed")

    nmatch <- namesmatchIts(obj1,obj2)

    if(matchnames && !nmatch) stop("names of the two inputs must match")
    if(overlap &&!overlapmatchesIts(obj1[,attr(nmatch,which="lut")],obj2)) stop("overlap data does not match")
    if(max(as.numeric(obj1@dates))<=max(as.numeric(obj2@dates))) {
        xlow <- obj1;    xhigh <- obj2
    } else {
        xlow <- obj2;    xhigh <- obj1
    }
    if(overlapsIts(obj1,obj2)) {
        highoverlap <- which(as.numeric(xhigh@dates)<=max(as.numeric(xlow@dates)) &
                             as.numeric(xhigh@dates)>=min(as.numeric(xlow@dates)))
        xhigh <- xhigh[-highoverlap,]
    }
    if(matchnames) {
        united <- rbind(xlow[,attr(nmatch,which="lut")],xhigh)
    } else {
        united <- rbind(xlow,xhigh)
    }

    alldates <- c(xlow@dates,xhigh@dates)
    result <- its(united,dates=alldates)
    return(result)
}
##daysecondIts-function----------------------------------------------
daysecondIts <- function (x,...)
{
    if (!inherits(x, "its")) stop("function is only valid for objects of class 'its'")
    hour <- as.POSIXlt(x@dates,...)$hour
    min <- as.POSIXlt(x@dates,...)$min
    sec <- as.POSIXlt(x@dates,...)$sec
    return(3600 * hour + 60 * min + sec)
}

##weekdayIts-function------------------------------------------------
weekdayIts <- function (x,...)
{
    if (!inherits(x, "its")) stop("function is only valid for objects of class 'its'")
    day <- as.POSIXlt(x@dates,...)$wday
    return((0 < day) & (day < 6))
}

##rangeIts-function--------------------------------------------------
rangeIts <- function(x,start=dates(x)[1],end=dates(x)[nrow(x)],format=its.format(),...)
{
    if (!inherits(x, "its")) stop("function is only valid for objects of class 'its'")
    if(mode(start)=="character") {start <- as.POSIXct(x=strptime(start,format=format),...)}
    if(mode(end)=="character") {end <- as.POSIXct(x=strptime(end,format=format),...)}
    if(!is(start,"POSIXt")) stop("start should be in POSIX format")
    if(!is(end,"POSIXt")) stop("end should be in POSIX format")
    start <- as.POSIXct(start)
    end <- as.POSIXct(end)
    return(x[which((x@dates>=start) & (x@dates<=end)),])
}

##newIts-function----------------------------------------------------
newIts <- function(x=NA,
                   start=format(Sys.Date(),format=its.format()),
                   end,
                   ncol=1,
                   by="DSTday",
                   extract=FALSE,
                   format=its.format(),
                   tz="",
                   ...)
{
    if(mode(start)=="character") {
        start.p <- as.POSIXct(x=strptime(start,format=format),tz=tz)
    } else {
        start.p <- as.POSIXct(start)
    }

    if(missing(end)) {
      end <- format(as.Date(start,format=its.format()) + 99,format=its.format())
      end.p <- as.POSIXct(x=strptime(end,format=format),tz=tz)
    } else if(mode(end)=="character") {
      end.p <- as.POSIXct(x=strptime(end,format=format),tz=tz)
    } else {
      end.p <- as.POSIXct(end)
    }

    dates <- seq(from=start.p,by=by,to=end.p)

    if(extract) {
        dates <- extractDates(dates=dates,...)
    }

    result <- its(matrix(x,ncol=ncol,nrow=length(dates)),dates)

    return(result)
}
##extractIts-function-----------------------------------------------
extractIts <- function(x,
                       weekday=FALSE,
                       find=c("all","last","first"),
                       period=c("week","month","year"),
                       partials=TRUE,
                       firstlast=FALSE,
                       select)
{
    if (!inherits(x, "its")) stop("function is only valid for objects of class 'its'")

    xdates <- extractDates(x@dates,
                           weekday=weekday,
                           find=find,
                           period=period,
                           partials=partials,
                           firstlast=firstlast,
                           select=select)

    return(x[dates=xdates])
}

##collapse-function---------------------------------------------------
collapseIts <- function(x)
{
    if (!inherits(x, "its")) stop("function is only valid for objects of class 'its'")

    labels <- dimnames(x)[[2]]
    uniquelabels <- unique(labels)
    jresult <- match(uniquelabels,labels)
    result <- x[,jresult]*NA

    if(length(uniquelabels)<length(labels)) {
        for(j in 1:length(uniquelabels)) {
            jall <- which(!is.na(match(labels,uniquelabels[j])))
            delta <- apply(x[,jall],1,var,na.rm=TRUE)
            if(any(delta[!is.na(delta)]>0))
                stop("column data must match in collapse function")

            result[,j] <- apply(x[,jall],1,median,na.rm=TRUE)
        }
    } else {
        result <- x
    }
    return(result)
}

##-Utilities-##############################################################################################################
##-Utility Methods-
##validity check----------------------------------------------------
validIts <-  function(object)
{
    if(!identical(nrow(object@.Data),length(object@dates)))
        return("Inconsistent length of dates")

    if(any(is.na(object@dates)))
        return("Missing values in dates")

    d <- diff(object@dates)
    if(any(d<0))
        return("Dates must be non-decreasing")

    ##if(!is.null(attr(object@dates,"tzone"))) return("Timezone attribute not allowed in dates slot of class its")
    return(TRUE)
}
setValidity("its",validIts)

##[-method-----------------------------------------------------------
subsetIts <- function(x, i, j, ..., drop=FALSE) {
    if(match("dates",names(list(...)),0)>0) {
        dates <- list(...)[["dates"]]
        if(!missing(i)) stop("cannot specify both dates and i")
        if(!is(dates,"POSIXt")) stop("dates should be in POSIX format")
        dates <- as.POSIXct(dates)
        i <- match(dates, dates(x))
        if(any(is.na(i))) stop("some dates are not found")
    }

    if(missing(drop)) {
        drop <- FALSE
    }

    if(missing(i)) {
        i <- min(1,nrow(x)):nrow(x)
    }

    if(missing(j)) {
        j <- min(1,ncol(x)):ncol(x)
    }

    subx <- x@.Data[i, j, drop = drop]
    dates <- x@dates[i]

    ans <- new("its",
               subx,
               dates=dates)

    return(ans)
}
setMethod("[", c("its","ANY"),subsetIts)

##"[<-"-method-----------------------------------------------------------
subset.replaceIts <- function(x, i, j, ..., value) {
    if(match("dates",names(list(...)),0)>0)
    {
        dates <- list(...)[["dates"]]
        if(!missing(i)) stop("cannot specify both dates and i")
        if(!is(dates,"POSIXt")) stop("dates should be in POSIX format")
        dates <- as.POSIXct(dates)
        i <- match(dates, dates(x))
        if(any(is.na(i))) stop("some dates are not found")
    }

    if(missing(i)) {
        i <- min(1,nrow(x)):nrow(x)
    }

    if(missing(j)) {
        j <- min(1,ncol(x)):ncol(x)
    }

    x@.Data[i, j] <- value@.Data
    x@dates[i] <- value@dates

    ans <- new("its",
               core(x),
               dates=dates(x))

    return(ans)
}
setReplaceMethod("[", signature(x="its", value="its"),subset.replaceIts)

##-Utility Functions-
##addDimnames-function-----------------------------------------------
addDimnames <- function(x) {
    if(is.null(dimnames(x))) {dimnames(x) <- list(NULL,NULL)}
    if(is.null(dimnames(x)[[1]])&(nrow(x)>0)) {dimnames(x)[[1]] <- 1:nrow(x)}
    if(is.null(dimnames(x)[[2]])&(ncol(x)>0)) {dimnames(x)[[2]] <- 1:ncol(x)}
    return(x)
}

##gapIts-function---------------------------------------------------
gapIts <- function(x,y,maxgap) {
    if (!inherits(x, "its")&inherits(y, "its"))
        stop("function is only valid for objects of class 'its'")

    d1 <- as.numeric(x@dates)
    d2 <- as.numeric(y@dates)
    dd1 <- diff(d1)
    dd2 <- diff(d2)
    del1 <- min(d1)-max(d2)
    del2 <- min(d2)-max(d1)
    gap <- max(max(del1,del2),0)
    allrange <- range(range(dd1),range(dd2))
    gapwrong <- gap>max(allrange[2],maxgap)
    attr(gap,which="gapexcessive") <- gapwrong
    return(gap)
}

##overlapsIts-function----------------------------------------------
overlapsIts <- function(x,y) {
    if (!inherits(x, "its") & inherits(y, "its"))
        stop("function is only valid for objects of class 'its'")

    if(max(as.numeric(x@dates))<max(as.numeric(y@dates))) {
        xlow <- x;
        xhigh <- y
    } else {
        xlow <- y;
        xhigh <- x
    }
    nooverlap <- max(as.numeric(xlow@dates))<min(as.numeric(xhigh@dates))
    return(!nooverlap)
}

##overlapmatchesIts-function----------------------------------------
overlapmatchesIts <- function(x,y) {
    if (!inherits(x, "its")&inherits(y, "its"))
        stop("function is only valid for objects of class 'its'")

    if(!overlapsIts(x,y)) {
        stop("no overlap")
    }

    if(max(as.numeric(x@dates))<=max(as.numeric(y@dates))) {
        xlow <- x;
        xhigh <- y
    } else {
        xlow <- y;
        xhigh <- x
    }

    if(min(as.numeric(xlow@dates))>min(as.numeric(xhigh@dates)))
        stop("appendor data must extend appendee data")

    lowoverlap <- which(as.numeric(xlow@dates)>=min(as.numeric(xhigh@dates)))
    highoverlap <- which(as.numeric(xhigh@dates)<=max(as.numeric(xlow@dates)))
    if(!identical(xlow@dates[lowoverlap],xhigh@dates[highoverlap])) {
        mymatch <- FALSE
    } else {
        mymatch <- identical(all.equal(xlow[lowoverlap,],xhigh[highoverlap,]),TRUE)
    }
    return(mymatch)
}

##namesmatchIts-function--------------------------------------------
namesmatchIts <- function(x,y) {
    lut <- match(dimnames(y)[[2]],dimnames(x)[[2]])
    namesmatch <- identical(dimnames(x[,lut])[[2]],dimnames(y)[[2]])&
    all(!duplicated(dimnames(x)[[2]]))&all(!duplicated(dimnames(y)[[2]]))
    attr(namesmatch,which="lut") <- lut
    return(namesmatch)
}

##its.format-function------------------------------------------------
its.format <- function(formatDefault=NULL) {
    if(is.null(formatDefault)) {
        outformat <- get(x=".itsformat",env=itsState,inherits=FALSE)
    } else {
        outformat <- formatDefault
        assign(x=".itsformat",value=formatDefault,env=itsState,inherits=FALSE)
    }
    return(outformat)
}
##expandIts-function------------------------------------------------
expandIts <- function(x)
{
    ##takes a single column 'its', if there are NAs, splits it into columns
    ##each column having only a single run of non-NA data
    if(all(is.na(x)))
        return(x)

    mat <- rbind(x[1,1,drop=FALSE]*NA,x[,1,drop=FALSE],x[nrow(x),1,drop=FALSE]*NA)
    ib <- which(diff(is.na(mat))==-1)
    ie <- which(diff(is.na(mat))==1)
    nruns <- length(ib)
    matexp <- matrix(NA,nrow=nrow(x),ncol=nruns)
    for(i in 1:nruns) {
        irun <- ib[i]:(ie[i]-1)
        matexp[irun,i] <- x[irun,1,drop=FALSE]
    }

    dimnames(matexp) <- list(dimnames(x)[[1]],rep(dimnames(x)[[2]][1],nruns))

    result <- its(matexp,x@dates)

    result
}

##extractDates-function----------------------------------------------
extractDates <- function(dates,
                         weekday=FALSE,
                         find=c("all","last","first"),
                         period=c("week","month","year"),
                         partials=TRUE,
                         firstlast=FALSE,
                         select)
{
    find <- match.arg(find)
    period <- match.arg(period)
    myindex1 <- 1:length(dates)

    ##1  optionally point only to weekdays
    if(weekday) {
        wday <- as.POSIXlt(dates)$wday
        myindex1 <- which( (0 < wday) & (wday < 6) )
    }

    if(period=="month") {
        theperiod <- 100*as.POSIXlt(dates[myindex1])$year+as.POSIXlt(dates[myindex1])$mon
        dayinperiod <- as.POSIXlt(dates[myindex1])$mday
    } else if(period=="year") {
        theperiod <- as.POSIXlt(dates[myindex1])$year
        dayinperiod <- as.POSIXlt(dates[myindex1])$yday
    } else if(period=="week") {
        theweek <- as.numeric(format(as.POSIXct(dates[myindex1]), "%U") )
        theyear <- as.numeric(format(dates[myindex1],"%Y"))
        incorrectPartialWeek <- theweek==0
        theyear[incorrectPartialWeek] <- theyear[incorrectPartialWeek]-1    ##first partial week in January assigned to last year
        theweek[incorrectPartialWeek] <- as.numeric(format(ISOdate(theyear[incorrectPartialWeek]-1,12,31),"%U"))    ##only incomplete Jan weeks are indexed 0 (see Jan 1995)
        theperiod <- 100*theyear+theweek
        dayinperiod <- as.POSIXlt(dates[myindex1])$wday
    }

    ##2  if selecting based on 'find'
    if(find=="all") {
        myindex2 <- 1:length(myindex1)
    } else {
        myindex2 <- setdiff(which(diff(c(theperiod[1],theperiod))!=0),1)
        if(find=="last") {
            myindex2 <- myindex2-1
        }
        if(partials) {
            if(find=="last") {
                myindex2 <- unique(c(myindex2,length(myindex1)))
            } else {
                myindex2 <- unique(c(1,myindex2))
            }
        }
    }

    ##3 select based on 'select'
    if(missing(select)) {
        myindex3 <- 1:length(myindex2)
    } else {
        myindex3 <- which(dayinperiod[myindex2]%in%select)
    }

    myindex <- myindex1[myindex2][myindex3]
    if(firstlast) {
        myindex <- unique(c(1,myindex,myindex1[length(myindex1)]))
    }

    if(all(is.na(myindex)))
        myindex <- NULL

    return(dates[myindex])
}

##gencurves-function----------------------------------------------
gencurves <- function(x)
{
    curves <- vector("list",(ncol(x)))
    for (j in 1:ncol(x)) {
        curves[[j]] <- list(x=as.numeric(x@dates),y=as.numeric(x[,j]@.Data))
    }

    return(curves)
}

##most.recent-function--------------------------------------------
most.recent <- function(x)
{
    ## return a vector of indices of the most recent TRUE value (thanks to Tony Plate)
    if (!is.logical(x)) stop("x must be logical")
    x.pos <- which(x)
    if (length(x.pos)==0 || x.pos[1] != 1) x.pos <- c(1, x.pos)
    rep(x.pos, c(diff(x.pos), length(x) - x.pos[length(x.pos)] + 1))
}
##locf-function---------------------------------------------------
locf <- function(x) {
    if (!inherits(x, "its"))
        stop("function is only valid for objects of class 'its'")
    y <- x
    jna <- which(apply(is.na(x),2,any))

    for(j in jna) {
        y[,j] <- y[most.recent(!is.na(y[,j])),j]
    }

    dates(y) <- dates(x)

    return(y)
}

##priceIts-function-------------------------------------------------
priceIts <- function (instruments = "^gdaxi",
                      start, end,
                      quote = c("Open","High", "Low", "Close"),
                      provider = "yahoo",
                      method = "auto",
                      origin = "1899-12-30",
                      compression="d",
                      quiet=TRUE)

    ## added new argument, compression
    ## may be "d", "w" or "m", for daily weekly or monthly
    ## John Bollinger, 2004-10-20, www.BollingerBands.com, bbands@yahoo.com
{
    if (provider != "yahoo")
        stop("provider not implemented")

    allinstruments <- NULL

    if (missing(start))
        start <- "1991-01-02"
    if (missing(end))
        end <- format(Sys.time() - 86400, "%Y-%m-%d")

    provider <- match.arg(provider)
    start <- as.POSIXct(start, tz = "GMT")
    end <- as.POSIXct(end, tz = "GMT")

    for(i in 1:length(instruments)) {
        url <- paste("http://chart.yahoo.com/table.csv?s=", instruments[i],
                     format(start, paste("&a=", as.character(as.numeric(format(start, "%m")) - 1), "&b=%d&c=%Y", sep = "")),
                     format(end,paste("&d=", as.character(as.numeric(format(end,"%m")) - 1), "&e=%d&f=%Y", sep = "")),
                     "&g=",compression,
                     "&q=q&y=0&z=", instruments[i], "&x=.csv", sep = "")

        destfile <- tempfile()
        status <- download.file(url, destfile, method = method, quiet=quiet)
        if (status != 0) {
            unlink(destfile)
            stop(paste("download error, status", status))
        }

        nlines <- length(count.fields(destfile, sep = "\n"))

        if (nlines == 1) {
            unlink(destfile)
            stop(paste("No data available for", instruments[i]))
        }
        data <- read.csv(destfile)
        data <- data[nrow(data):1,]         ## and inverse order in data
        y <- its(as.matrix(data[,-1]),
                 dates=strptime(as.character(data[,1]), format="%Y-%m-%d"))

        oneinstrument <- its(y)[,quote]

        names(oneinstrument) <- paste(instruments[i],quote)
        allinstruments <- union(allinstruments,oneinstrument)
    }
    return(allinstruments)
}
