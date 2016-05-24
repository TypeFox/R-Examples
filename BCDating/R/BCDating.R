BCDating <- setClass("BCDating",
         representation(name="character", 
                        states="ts",
                        peaks="numeric",
                        troughs="numeric",
                        y="ts",
                        param="list",
                        type="character"))
# Functions on BBQ -------------------------

BBQ <- function (y,
                 mincycle = 5, minphase = 2, 
                 name = "") 
{
  k.peak = 2
  k.trough = 2
  l.peak = 2
  l.trough = 2
  e = NULL
  if (!is.ts(y)) 
    stop("Argument <y> should be an object of class ts")
  e1 <- BCDating.init(y, k.peak = k.peak, k.trough = k.trough, 
                      l.peak = l.peak, l.trough = l.trough)
  datok <- BCDating.censor(e1, y, mincycle = mincycle, 
                           minphase = minphase, e = e)
  datok@y <- y
  datok@type <- "BCDating:BBQ"
  datok@param <- list(k.peak = k.peak, k.trough = k.trough,
                      l.peak = l.peak, l.trough = l.trough, 
                      mincycle = mincycle, minphase = minphase, e = e)
  return(datok)
}

build.mat_tp <- function (peaks, troughs) 
{
    if (any(is.na(peaks), is.na(troughs))) 
        stop("Missing values are not allowed")
    np <- length(peaks)
    nt <- length(troughs)
    peaksn <- matrix(peaks, np, 2)
    peaksn[, 2] <- 1
    troughsn <- matrix(troughs, nt, 2)
    troughsn[, 2] <- 0
    anmat <- rbind(peaksn, troughsn)
    anmat <- anmat[order(anmat[, 1]), ]
    return(anmat)
}

CTS_BBQ <-function (x, i) 
return(ETS_BBQ(-x, i))

ETS_BBQ <- function (x, i) 
{
  if (!((i >= 1) & (i <= length(x)))) 
    stop("wrong parameter i")
  return(all(x[i] >= x))
}

BCDating.alter1 <- function (dat, y) 
{
    if (class(dat)[1] != "BCDating") 
        stop("Argument <dat> must be of an object of class 'BCDating'")
    bcp <- dat@peaks
    bct <- dat@troughs
    np <- length(bcp)
    nt <- length(bct)
    keep <- rep(TRUE, nt)
    if (bcp[1] < bct[1]) { # Peaks are first
        r <- min(np, nt)
        for (i in 1:r) 
          if (y[bct[i]] > y[bcp[i]]) 
            keep[i] <- FALSE
    }
    else { # Troughs are first
        r <- min(np, nt - 1)
        for (i in 1:r) if (y[bct[i + 1]] > y[bcp[i]]) 
            keep[i + 1] <- FALSE
    }
    if (r > nt) 
        keep <- keep[1:nt]
    bct <- bct[keep]
    res <- dat
    res@peaks <- bcp
    res@troughs <- bct
    return(res)
}

BCDating.alter2 <- function (dat, y)
{
    if (class(dat)[1] != "BCDating") 
        stop("Argument <dat> must be of an object of class 'BCDating'")
    bcp <- dat@peaks # Business Cycle Peaks
    bct <- dat@troughs # Business Cycle Throughs
    if (any(is.na(bcp), is.na(bct), is.na(y))) 
        stop("Missing values are not allowed")
    np <- length(bcp) # Number of Peaks
    nt <- length(bct) # Number of Throughs
    anmat <- build.mat_tp(bcp, bct)
    j <- 1
    repeat {
        state1 <- anmat[j, 2]
        state2 <- anmat[j + 1, 2]
        if (state1 == state2) {
            if (state1 == 1) {
                if (y[anmat[j, 1]] > y[anmat[j + 1, 1]]) 
                  vire <- j + 1
                else vire <- j
            }
            if (state1 == 0) {
                if (y[anmat[j, 1]] > y[anmat[j + 1, 1]]) 
                  vire <- j
                else vire <- j + 1
            }
            anmat <- anmat[-vire, ]
        }
        else j <- j + 1
        if (j >= nrow(anmat)) 
            break
    }
    if (FALSE) {
        res <- dat
        res@peaks <- anmat[anmat[, 2] == 1, 1]
        res@troughs <- anmat[anmat[, 2] == 0, 1]
        res@states <- BCDating.pt2states(start = start(y), end = end(y), 
            freq = frequency(y), peaks = res@peaks, troughs = res@troughs)
    }
    res <- BCDating.peakstroughs(start = start(y), end = end(y), 
        freq = frequency(y), peaks = anmat[anmat[, 2] == 1, 1], 
        troughs = anmat[anmat[, 2] == 0, 1], name = dat@name, 
        type = dat@type, param = dat@param)
    return(res)
}

BCDating.censor <- function (dat, y, mincycle = 5, minphase = 2, e = NULL) 
{
    if (class(dat)[1] != "BCDating") 
        stop("argument <dat> should be an object of class 'BCDating'")
    if (class(y)[1] != "ts") 
        stop("argument <y> should be an object of class ts'")
    if (!(frequency(y) %in% c(4, 12))) 
        stop("the time series (argument <y>) should be quarterly or monthly")
    if (is.null(e)) {
        if (frequency(y) == 12) 
            e <- 6
        if (frequency(y) == 4) 
            e <- 2
    }
    et1 <- BCDating.alter2(dat, y)
    et2 <- BCDating.alter1(et1, y)
    et3 <- BCDating.alter2(et2, y)
    deb <- et3
    repeat {
        init <- deb
        et1 <- BCDating.enf1p(init, y, mincycle = mincycle)
        et2 <- BCDating.enfvbp(et1, y, e = e)
        et3 <- BCDating.enfvc(et2, y)
        et4 <- BCDating.enf1p(et3, y, mincycle = mincycle)
        et5 <- BCDating.enfve(et4, y, minphase = minphase)
        et6 <- BCDating.enfvc(et5, y)
        fin <- et6
        if (all(deb@states == fin@states)) 
            break
        deb <- fin
    }
    return(fin)
}

BCDating.enf1p <- function (dat, y, mincycle = 5) 
{
    if (class(dat)[1] != "BCDating") 
        stop("argument <dat> must be of an object of class 'BCDating'")
    bcp <- dat@peaks
    i <- 2
    repeat {
        if (bcp[i] - bcp[i - 1] < mincycle) {
            if (y[bcp[i]] > y[bcp[i - 1]]) 
                vire <- i - 1
            else vire <- i
            bcp <- bcp[-vire]
        }
        else i <- i + 1
        if (i >= length(bcp)) 
            break
    }
    dat@peaks <- bcp
    intermediaire <- BCDating.alter2(dat, y)
    bct <- intermediaire@troughs
    i <- 2
    repeat {
        if (bct[i] - bct[i - 1] < mincycle) {
            if (y[bct[i]] < y[bct[i - 1]]) 
                vire <- i - 1
            else vire <- i
            bct <- bct[-vire]
        }
        else i <- i + 1
        if (i >= length(bct)) 
            break
    }
    dat@troughs <- bct
    fin <- BCDating.alter2(dat, y)
    return(fin)
}

BCDating.enfvbp <- function (dat, y, e = 6) 
{
    if (class(dat)[1] != "BCDating") 
        stop("argument <dat> should be an object of class 'BCDating'")
    if (class(y)[1] != "ts") 
        stop("argument <y> should be an object of class ts'")
    n <- length(y)
    proceed <- function(seqq, e) return(seqq[(seqq > e) & (seqq <= 
        n - e)])
    res <- dat
    res@peaks <- proceed(dat@peaks, e)
    res@troughs <- proceed(dat@troughs, e)
    return(BCDating.alter2(res, y))
}

BCDating.enfvc <- function (dat, y) 
{
    if (class(dat)[1] != "BCDating") 
        stop("argument <dat> should be an object of class 'BCDating'")
    if (class(y)[1] != "ts") 
        stop("argument <y> should be an object of class ts'")
    bcp <- dat@peaks
    bct <- dat@troughs
    n <- length(y)
    repeat {
        nothing_done <- TRUE
        if ((length(bcp) == 0) | (length(bct) == 0)) 
            break
        m <- min(min(bcp), min(bct))
        if (m == min(bcp)) {
            change_p <- TRUE
            change_t <- FALSE
            if (y[1] > y[bcp[1]]) 
                bcp <- bcp[-1]
            else change_p <- FALSE
        }
        else {
            change_p <- FALSE
            change_t <- TRUE
            if (y[1] < y[bct[1]]) 
                bct <- bct[-1]
            else change_t <- FALSE
        }
        nothing_done <- nothing_done & (!change_p) & (!change_t)
        m <- max(max(bcp), max(bct))
        if (m == max(bcp)) {
            np <- length(bcp)
            change_p <- TRUE
            change_t <- FALSE
            if (y[n] > y[bcp[np]]) 
                bcp <- bcp[-np]
            else change_p <- FALSE
        }
        else {
            nt <- length(bct)
            change_p <- FALSE
            change_t <- TRUE
            if (y[n] < y[bct[nt]]) 
                bct <- bct[-nt]
            else change_t <- FALSE
        }
        nothing_done <- nothing_done & (!change_p) & (!change_t)
        if (nothing_done) 
            break
    }
    res <- dat
    res@peaks <- bcp
    res@troughs <- bct
    return(BCDating.alter2(res, y))
}

BCDating.enfve <- function (dat, y, minphase = 2) 
{
    if (class(dat)[1] != "BCDating") 
        stop("argument <dat> must be of an object of class 'BCDating'")
    j <- 1
    repeat {
        anmat <- build.mat_tp(dat@peaks, dat@troughs)
        if (j >= nrow(anmat)) 
            break
        if ((anmat[j + 1, 1] - anmat[j, 1]) < minphase) {
            anmat <- anmat[-(j + 1), ]
            dat@peaks <- anmat[anmat[, 2] == 1, 1]
            dat@troughs <- anmat[anmat[, 2] == 0, 1]
            dat <- BCDating.alter2(dat, y)
        }
        else j <- j + 1
    }
    return(dat)
}

BCDating.init <- function (y, ETS = ETS_BBQ, CTS = CTS_BBQ, 
                           k.peak = 2, k.trough = 2, 
                           l.peak = 2, l.trough = 2) 
{
    if (!(is.ts(y))) 
        stop("Argument <y> should be an object of class 'ts'")
    n <- length(y)
    peaks <- rep(NA, n)
    troughs <- rep(NA, n)
    for (i in 1:n) {
        LB_p <- max(1, i - k.peak)
        LB_t <- max(1, i - k.trough)
        z_p <- y[LB_p:min(n, i + l.peak)]
        z_t <- y[LB_t:min(n, i + l.trough)]
        if ((i > k.peak) & (i <= n - l.peak)) 
            peaks[i] <- ETS(z_p, i - LB_p + 1)
        if ((i > k.trough) & (i <= n - l.trough)) 
            troughs[i] <- CTS(z_t, i - LB_t + 1)
    }
    return(BCDating.peakstroughs(start = start(y), end = end(y),
                                 freq = frequency(y), 
                                 peaks = which(peaks), 
                                 troughs = which(troughs)))
}

BCDating.peakstroughs <- function (start, end, freq = NULL, 
                                   peaks, troughs, name = "", 
                                   type = "user-defined", param = NULL) 
{
    if (is.null(freq)) {
        freq <- 0
        if (substr(peaks[1], 5, 5) == "M") 
            freq <- 12
        if (substr(peaks[1], 5, 5) == "Q") 
            freq <- 4
        char2time <- function(chaine) {
            year <- as.integer(substr(chaine, 1, 4))
            per <- as.integer(substr(chaine, 6, nchar(chaine)))
            return(year + (per - 1)/freq)
        }
        temps <- ts(0, start = start, end = end, frequency= freq)
        temps <- time(temps)
        peaks <- sapply(char2time(peaks), function(a) which(abs(a - 
            temps) < 0.001))
        troughs <- sapply(char2time(troughs), function(a) which(abs(a - 
            temps) < 0.001))
    }
    else {
        temps <- ts(0, start = start, end = end, frequency = freq)
        temps <- time(temps)
    }
    if (!(freq %in% c(4, 12))) 
        stop("frequency must be 12 (monthly dates) or 4 (quarterly dates)")
    states <- BCDating.pt2states(start, end, freq, peaks, troughs)
    if (is.null(param)) 
        param <- vector("list", 0)
    return((new("BCDating", name = name, states = states, peaks = peaks, 
        troughs = troughs, param = param, type = type)))
}

BCDating.pt2states <- function (start, end, freq, peaks, troughs) 
{
    states <- ts(0, start = start, end = end, frequency = freq)
    n <- length(states)
    mat_tp <- build.mat_tp(peaks, troughs)
    r <- nrow(mat_tp)
    if (mat_tp[r, 1] < n) 
        mat_tp <- rbind(mat_tp, c(n, 1 - mat_tp[r, 2]))
    if (peaks[1] < troughs[1]) 
        add <- 0
    else add <- 1
    states[1:mat_tp[1, 1]] <- (-1)^add
    for (j in 1:(nrow(mat_tp) - 1)) 
    {
      states[(mat_tp[j, 1] + 1):mat_tp[j + 1, 1]] <- (-1)^(j + add)
    }
    return(states)
}

# Show Method ----------------------------

matsummary <- function (object) 
{
    np <- length(object@peaks)
    nt <- length(object@troughs)
    r <- max(np, nt)
    if (r != 0) {
        res <- matrix(NA, r, 2)
        if (np == 0) 
            res[1, 2] <- object@troughs[1]
        if (nt == 0) 
            res[1, 1] <- object@peaks[1]
        if ((np > 0) & (nt > 0)) {
            if (object@peaks[1] < object@troughs[1]) {
                res[1:np, 1] <- object@peaks
                res[1:nt, 2] <- object@troughs
            }
            else {
                if (np == nt) 
                  res <- matrix(NA, r + 1, 2)
                res[2:(np + 1), 1] <- object@peaks
                res[1:nt, 2] <- object@troughs
            }
        }
    }
    colnames(res) <- c("Peaks", "Troughs")
    return(res)
}

matsummary2 <- function (object) 
{
    dat <- object
    nr <- length(dat@peaks) + length(dat@troughs) + 1
    summarytab <- matrix(NA, nr, 7)
    change.states <- sort(c(0, dat@peaks, dat@troughs))
    summarytab[,1] <- dat@states[change.states + 1]
    summarytab[,3] <- c(change.states[-1], NA)
    summarytab[-1,2] <- summarytab[-nr,3]
    summarytab[,4] <- summarytab[,3] - summarytab[,2]
    summarytab[,5] <- dat@y[summarytab[,2]]
    summarytab[,6] <- dat@y[summarytab[,3]]
    summarytab[,7] <- summarytab[,1] * (summarytab[,6] - summarytab[,5])
    return(summarytab)
}

if (!isGeneric("show")) {
    setGeneric("show", function(object, ...) standardGeneric("show"))
 }

setMethod("show",
    signature(object = "BCDating"),
    function (object) 
    {
      if (nchar(object@name) > 0) 
        cat("Dating name :", object@name, "\n")
      res <- matsummary(object)
      affich <- res
      affich[, c(1, 2)] <- ts2char(object@states)[res]
      duration <- res[, 2] - res[, 1]
      affich <- cbind(affich, duration)
      colnames(affich)[3] <- "Duration"
      res <- data.frame(affich)
      print(res)
    }
          )

# Plot Method ---------------------------------------------

vrt <- function(v) #Virtual Time series to prepare the plot
{
  st <- min(time(v),na.rm=TRUE)
  end <- max(time(v),na.rm=TRUE)
  dl <- length(v)
  tsmin <- as.numeric(min(v,na.rm=TRUE))
  tsmax <- as.numeric(max(v,na.rm=TRUE))
  tsr <- tsmax-tsmin
  vrt <- ts(data=1:dl/dl*tsr*1.1+tsmin-0.05*tsr,start=st,frequency=frequency(v))
}

if (!isGeneric("plot")) {
    setGeneric("plot", function(x,y, ...) standardGeneric("plot"))
 }

setMethod("plot",
    signature(x = "BCDating", y = "missing"),
    function (x, y, dates =  FALSE, yearrep=2,
              col.bg=grey(0.8),col.exp=grey(1),col.rec = grey(0.45),
              xaxs="i", yaxs="i",
              main = "", xlab = "", ylab = "", lwd=2, cex = 0.5, 
              vert=NULL, col.vert="darkblue",
              xmin=NULL, xmax=NULL, ymin=0, ymax=1,
              debug=FALSE, ...) 
    {
      if(!is.null(main))
        if ((nchar(x@name) > 0) & (nchar(main) == 0)) 
          main <- paste(x@name)
      
      smin <- as.numeric(min(time(x@states),na.rm=TRUE))
      smax <- as.numeric(max(time(x@states),na.rm=TRUE))
      f <- frequency(x@states)
      
      if(is.null(xmin))
      {
        xmin <- smin - 1/f
        xmin <- round(xmin*f)/f
      }
      if(is.null(xmax))
      {      
        xmax <- smax + 1/f
        xmax <- round(xmax*f)/f
      }
             
    suppressWarnings(
      xstates <- window(ts(c(rep(NA,f*(smin-xmin)),x@states,rep(NA,f*(xmax-smax))),
                        start=xmin,frequency=f),xmin,xmax))
    xpeaks <- x@peaks + (smin-xmin)*f
    xtroughs <- x@troughs + (smin-xmin)*f
    
#    xn <- length(xstates)


    plot(c(xmin,xmax),c(ymin,ymax),type="n", xaxs = xaxs, yaxs = yaxs, 
         main = "", xlab = "", ylab = "",...)
    
    #   back.color <- ts(data=rep(col.bg, xn),start=xmin,frequency=f) #no info
    # 
    #   back.color[which(xstates == -1)] <- col.rec #recession
    #   back.color[which(xstates != -1)] <- col.exp #expansion
    
    temps <- time(xstates)
    dt <- deltat(xstates)
    #   rect(temps - 0.5 * dt, ymin, temps + 0.5 * dt, ymax, col = back.color,border = NA)
    
    ## Plotting Recessions
    if(xpeaks[1]>xtroughs[1]){
      a <- (smin-xmin)*f+1
    }else{
      a <- NULL
    }    
    if(xpeaks[length(xpeaks)]>xtroughs[length(xtroughs)]){
      b <- (smax-xmin)*f+1
    }else{
      b <- NULL
    }
    rbxs <- temps[c(a,xpeaks+1)]
    rbxe <- temps[c(xtroughs,b)]
    
    rect(rbxs-0.5*dt,ymin,rbxe+0.5*dt,ymax,col=col.rec,border=NA)
    
    ## Plotting Expansions
    if(xpeaks[1]<xtroughs[1]){
      a <- (smin-xmin)*f+1
    }else{
      a <- NULL
    }    
    if(xpeaks[length(xpeaks)]<xtroughs[length(xtroughs)]){
      b <- (smax-xmin)*f+1
    }else{
      b <- NULL
    }
    
    ebxs <- temps[c(a,xtroughs+1)]
    ebxe <- temps[c(xpeaks,b)]
    
    rect(ebxs-0.5*dt,ymin,ebxe+0.5*dt,ymax,col=col.exp,border=NA)

    ## Plotting Unknowns
    ubxs <- c(temps[1],smax+1/f)
    ubxe <- c(smin-1/f,temps[length(temps)])
    
    rect(ubxs-0.5*dt,ymin,ubxe+0.5*dt,ymax,col=col.bg,border=NA)
    
    title(main=main,ylab=ylab,xlab=xlab)
    if(dates)
      {
        Dates <- paste(substr(time(xstates),5-yearrep,4)
                     ,cycle(xstates),sep=":")
#        show(Dates) # ############ JUST FOR DEBUG
       for(p in xpeaks)
       {
         text(x=time(xstates)[p]+.5/f,y=0.8*ymax,
             labels=Dates[p],pos=4,offset=0,cex=cex)
         lines(x=rep(time(xstates)[p],2),y=c(0.7,0.9),col="blue")
        }
      for(p in xtroughs)
       {
         text(x=time(xstates)[p]+.5/f,y=0.2*ymax,
             labels=Dates[p],pos=4,offset=0,cex=cex)
       }
     }
    if(!is.null(vert))
      segments(x0=vert,y0=0,x1=vert,y1=1,col=col.vert,lwd=3,lty=2)
    box(which = "plot")
    }
)

setMethod("plot",
          signature(x = "BCDating", y = "ts"),
          function (x, y,main = "",
                    window=FALSE,Dwindow=FALSE,averages=FALSE,dates=FALSE,yearrep=2,
                    col="red",col.bg=grey(.8),col.exp=grey(1),col.rec=grey(.45),
                    cex = 0.5,xlab = "", ylab = "",lwd=2, 
                    vert=NULL, col.vert="darkblue",
                    xmin=NULL, xmax=NULL, ymin=0, ymax=1,
                    debug=FALSE, ...) 
          {
            if(averages & !window)
            {
              warning("BCDating: Plotting Averages only in windowed mode")
              window = TRUE
            }
            
            smin <- as.numeric(min(time(x@states),na.rm=TRUE))
            smax <- as.numeric(max(time(x@states),na.rm=TRUE))
            f <- frequency(x@states)
            
            if(is.null(xmin))
            {
              xmin <- smin - 1/f
              xmin <- round(xmin*f)/f
            }
            if(is.null(xmax))
            {      
              xmax <- smax + 1/f
              xmax <- round(xmax*f)/f
            }
            
            suppressWarnings(
              xstates <- window(ts(c(rep(NA,f*(smin-xmin)),x@states,rep(NA,f*(xmax-smax))),
                                   start=xmin,frequency=f),xmin,xmax))
            xpeaks <- x@peaks + (smin-xmin)*f
            xtroughs <- x@troughs + (smin-xmin)*f
            
            if(window & !Dwindow)
            {
              ymin <- as.numeric(min(window(y,start=tsp(x@states)[1],end=tsp(x@states)[2]),na.rm=TRUE))
              ymax <- as.numeric(max(window(y,start=tsp(x@states)[1],end=tsp(x@states)[2]),na.rm=TRUE))
            }else if(!window & !Dwindow)
            {
              xmin <- as.numeric(min(time(xstates),time(y),na.rm=TRUE))
              xmax <- as.numeric(max(time(xstates),time(y),na.rm=TRUE))
              ymin <- as.numeric(min(y,na.rm=TRUE))
              ymax <- as.numeric(max(y,na.rm=TRUE))
            }else if(window & Dwindow)
            {
              xminy <- min(time(y),na.rm=TRUE)
              xmind <- min(time(xstates),na.rm=TRUE)
              xmin <- as.numeric(max(xminy,xmind,na.rm=TRUE))
              xmaxy <- max(time(y),na.rm=TRUE)
              xmaxd <- max(time(xstates),na.rm)
              xmax <- as.numeric(min(xmaxy,xmaxd),na.rm=TRUE)
              ymin <- as.numeric(min(window(y,start=xmin,end=xmax),na.rm=TRUE))
              ymax <- as.numeric(max(window(y,start=xmin,end=xmax),na.rm=TRUE))
            }else if(!window & Dwindow)
            {
              xmin <- as.numeric(min(time(y),na.rm=TRUE))
              xmax <- as.numeric(max(time(y),na.rm=TRUE))
              ymin <- as.numeric(min(y,na.rm=TRUE))
              ymax <- as.numeric(max(y,na.rm=TRUE))
            }
            
            suppressWarnings(
              xstates <- window(ts(c(rep(NA,f*(smin-xmin)),x@states,rep(NA,f*(xmax-smax))),
                                   start=xmin,frequency=f),xmin,xmax))
              
            yl <- ymax-ymin
            ymin <- ymin - yl/40
            ymax <- ymax + yl/40
            plot(x,dates=dates,yearrep=yearrep,
                 col.bg=col.bg,col.exp=col.exp, col.rec=col.rec,
                 main=main, xlab=xlab, ylab=ylab, lwd=lwd, cex=cex,
                 xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax,
                 vert=vert, col.vert=col.vert, debug=debug, ...)
            y <- cbind(y,0)
            if(length(averages)==1) 
              averages <- rep(averages,dim(y)[2]-1)
            if(length(col)==1 & !is.null(dim(y)))
              col <- rep(col,dim(y)[2]-1)
            
            for(v in 1:(dim(y)[2]-1)){
              if(averages[v]==FALSE)
                lines(y[,v],lwd=lwd,col = col[v])
              else{
                a <- avgts(y[,v],x)
                l <- length(x@states)
                pt <- c(1,x@troughs,x@peaks,l)
                pt <- pt[order(pt)]
                lpt <- length(pt)
                spt <- pt[1:lpt-1]
                ept <- pt[2:lpt]
                add <- 0.5/frequency(x@states)
                segments(time(a)[spt]-add,a[spt+1],time(a)[ept]-add,a[ept-1],lwd=lwd,col=col[v])
              }
            }
            #   points(y, pch = pch, col = col[2], cex = cex)
            if(!is.null(vert))
              segments(x0=vert,y0=ymin,x1=vert,y1=ymax,col=col.vert,lwd=3,lty=2)
            
            box(which = "plot")
          }
)


setMethod("plot",
    signature(x = "ts", y = "BCDating"),
    function (x, y, ...) 
    {
        plot(y,x,...)
    }
)

setMethod("plot",
    signature(x = "BCDating", y = "BCDating"),
    function (x, y, ...) 
    {
        plot(list(x,y),...)
    }
)

setMethod("plot",
          signature(x = "list", y = "missing"),
          function (x,pch = 1,cex = 0.8,dates=TRUE,yearrep=2,lines=4, ...) 
{
  if (class(x[[1]])[1] != "BCDating")
    stop("argument <x> should be of an array of objects of class 'BCDating'")
  xax <- x[[1]]@states
  rx <- range(time(xax))
  rx <- c(rx[1],rx[2]+0.5)
  ry <- c(0,length(x))
  opar <- par(mar=c(lines,4,4,2)+0.1)
  plot(rx, ry, type = "n", xaxs = "i", yaxt = "n", 
       xlab = "", ylab = "", ...)
  par(opar)
  
  for(i in 1:length(x))
  {
    if (class(x[[i]])[1] != "BCDating") 
      stop("argument <x> should be of an array of objects of class 'BCDating'")
    back.color <- rep(NA, length(x[[i]]@states))
    back.color[x[[i]]@states == -1] <- grey(0.4+0.3*i/length(x))
    back.color[x[[i]]@states != -1] <- grey(1)
    temps <- time(x[[i]]@states)
    dt <- deltat(x[[i]]@states)
    rect(temps , i-1, temps + dt, i, col = back.color, 
         border = NA)
    
    if(dates)
    {
      Dates <- paste(substr(time(x[[i]]@states),5-yearrep,4)
                     ,cycle(x[[i]]@states),sep=":")
      for(p in x[[i]]@peaks)
      {
        text(x=time(x[[i]]@states)[p]-0.1,y=i-0.2,
             labels=Dates[p],pos=4,cex=cex)
      }
      for(p in x[[i]]@troughs)
      {
        text(x=time(x[[i]]@states)[p]-0.1,y=i-0.8,
             labels=Dates[p],pos=4,cex=cex)
      }
    }
  }
}
)

# Summary Method --------------------------------

if (!isGeneric("summary")) {
    setGeneric("summary", function(object, ...) standardGeneric("summary"))
 }

setMethod("summary",
    signature(object = "BCDating"),
    function (object, print=TRUE, ...) 
    {
      summarytab <- matsummary2(object)
      indic <- matrix(NA, 2, 2)
      colnames(indic) <- c("Amplitude", "Duration")
      rownames(indic) <- c("Exp=]T;P]", "Rec=]P;T]")
      indic[1, 1] <- mean(summarytab[summarytab[, 1] == 1, 7],na.rm = TRUE)
      indic[2, 1] <- mean(summarytab[summarytab[, 1] == -1, 7],na.rm = TRUE)
      indic[1, 2] <- mean(summarytab[summarytab[, 1] == 1, 4],na.rm = TRUE)
      indic[2, 2] <- mean(summarytab[summarytab[, 1] == -1, 4],na.rm = TRUE)
      if (isTRUE(print)) {
          df.print <- as.data.frame(summarytab)
          df.print[summarytab[, 1] == 1, 1] <- "Expansion"
          df.print[summarytab[, 1] == -1, 1] <- "Recession"
          df.print[, c(2, 3)] <- ts2char(object@states)[summarytab[, 
              c(2, 3)]]
          df.print[, c(5, 6)] <- round(summarytab[, c(5, 6)], 0)
          df.print[, 7] <- round(summarytab[, 7], 1)
          colnames(df.print) <- c("Phase", "]Start", ";End]", "Duration", 
              "LevStart", "LevEnd", "Amplitude")
          print(df.print)
          cat("\n")
          print(round(indic, 1))
      }
      return(invisible(indic))
    }
)

ts2char <- function (obj) 
{
    if (!is.ts(obj)) 
        stop("argument <obj> should be an object of class 'ts'")
    temps <- time(obj)
    year <- floor(temps)
    per <- cycle(obj)
    letter <- "."
    if (frequency(obj) == 4) 
        letter <- "Q"
    if (frequency(obj) == 12) 
        letter <- "M"
    return(paste(year, letter, per, sep = ""))
}

avgmat <- function (Dating) {
  
  l <- length(Dating@states)
  points <- c(Dating@troughs,Dating@peaks)
  points <- points[order(points)]
  x <- matrix(rep.int(0,l*l),l,l)
  x[1:points[1],1:points[1]] <- 1/points[1]
  for(i in 1:(length(points)-1))
    x[(points[i]+1):points[i+1],(points[i]+1):points[i+1]] <- 1/(points[i+1]-points[i])
  if(points[length(points)]+1 <= l)
    x[(points[length(points)]+1):l,(points[length(points)]+1):l] <- 1/(l-points[length(points)])
  avgmat <- x
}

avgts <- function(ts,Dating)
{
  ts <- na.omit(ts)
  wts <- window(ts,start=start(Dating@states),end=end(Dating@states),
                frequency=frequency(Dating@states))
  wdat <- window(Dating,start=start(wts),end=end(wts))
  avgts <- ts(avgmat(wdat) %*% wts,start=start(wts),
              frequency=frequency(wts))
}
          
sid <- function(ts1,ts2)
{
  sd <- start(ts2)-start(ts1)
  sid <- sd[1]*frequency(ts1)+sd[2]
}

# window Method ----------------------------------------------------

if (!isGeneric("window")) {
    setGeneric("window", function(x, ...) standardGeneric("window"))
 }
          
setMethod("window",
    signature(x = "BCDating"),
    function (x, ...)
    {
      d <- x
      d@states <- window(x@states,...)
      d@peaks <- x@peaks + sid(d@states,x@states)      
      d@troughs <- x@troughs + sid(d@states,x@states)
      d@peaks <- d@peaks[d@peaks <= length(d@states) & d@peaks >=1]
      d@troughs <- d@troughs[d@troughs <= length(d@states) & d@troughs >=1]
      d@y <- window(x@y,...)
      return(d)
    }
)
