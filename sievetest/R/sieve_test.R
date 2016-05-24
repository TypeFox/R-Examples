#
#  Copyright 2014 Petr Matoušů <pmatousu@more-praha.cz>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.


#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#:  desc.std    makes description metadata
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

desc.std <- function(Description=NA,Title=NA,Sample=NA,Unit=NA,Date=NA,x=NULL)  {

#  merges lists, if not NA y overwrites x
ml <- function(x, y) {
  x.names <- names(x)
  y.names <- names(y)
  m.names <- sort(unique(c(x.names, y.names)))
  sapply(m.names, function(i) {
    if (is.list(x[[i]]) & is.list(y[[i]])) ml(x[[i]], y[[i]])
    else if (i %in% y.names & !is.na(y[[i]])) y[[i]]
    else x[[i]]
  }, simplify = FALSE)
}

  new <- list(Description=Description,
              Title=Title,
              Sample=Sample,
              Unit=Unit,
              Date=Date)

  if(is.null(x)) {
    ret <- new
  } else {
    stopifnot(is.std(x),length(x) == 1)
    old <- list(Description=x[[1]]$Description,
                Title=x[[1]]$Title,
                Sample=x[[1]]$Sample,
                Unit=x[[1]]$Unit,
                Date=x[[1]]$Date)
    ret <- ml(old,new)
  }

  return(lapply(ret,as.character))
}


#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#:  create std object by hand
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

std <- function(a,r,desc=desc.std(),lmargs=list())  {
  # a ... apertures
  # r ... mass % retained on sieve
  stopifnot(length(a) == length(r),is.list(lmargs),is.list(desc))
  tab <- data.frame(a,r)
  tab <- na.exclude(tab)
  metadata <- desc
  names(tab) <- c("a","r")
  if (!sum(tab$r) == 100) {
    warning("Sum of retained is not 100, recalculating to percents, check results.")
    tab$r <- 100*tab$r/sum(tab$r)
  }
  if (!all(tab$r > 0)) {
    warning("Zero entries in retained field, excluding, check results.")
    ind <- tab$r > 0
    tab <- tab[ind,]
  }
  if (all(sapply(desc[-1],is.na)))
    warning("Insufficient metadata.")

  # Ordering
  tab <- tab[order(tab$a,decreasing=T),]

  # Cumulative Oversize
  tab$O <- cumsum(tab$r)
  # Based on approach, the last number of tab$O should be 100,
  # as WAR of R acurracy, the 100 is forced
  # to avoid log(log(100/100.00000000000001))  problem (NaNs produced)
  tmp <- rev(tab$O)
  if(tmp[1] - 100 < 1e-12) {
    tmp[1] <- 100
    tab$O <- rev(tmp)
  }

  # Undersize
  tab$U <- 100 - tab$O
  tmp <- rev(tab$U)
  if(tmp[1] < 1e-12) {
    tmp[1] <- 0
    tab$U <- rev(tmp)
  }

  # log x
  tab$X <- x2rrx(tab$a)
  # log log y
  tab$Y <- y2rry(tab$O)
  # Linear fit
  l <- do.call('lm', c(list(formula=Y ~ X, data=tab[tab$a > 0,]),lmargs))
  # RR coefficients
  xs <- rrx2x((-1)*l$coefficients[1]/l$coefficients[2])
  names(xs) <- "xs"
  ex <- l$coefficients[2]
  names(ex) <- "ex"
  RRcoeff <- list(xs=xs,ex=ex)
  names(RRcoeff) <- c("diameter","exponent")
  # size features
  a_modus <- unname(xs*((ex-1)/ex)^(1/ex))
  a_mean <- unname(xs*(log(2))^(1/ex))
  o90um <- unname(orr(90,ex,xs))
  size <- list(modus=a_modus,mean=a_mean,o90um=o90um)
  ret <- list(stdata=tab,lmfit=l,RRcoefficients=RRcoeff,size=size)
  ret <- list(c(metadata,ret))

  # class
  class(ret) <- c("std",class(ret))
  return(ret)
}


#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#:  read std data from CSV
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

read.std <- function(file, sep="\t", dec=".")  {
  stopifnot(file.exists(file))
  # get number of columns
  hd <- ncol(read.table(file,header=F,sep=sep,stringsAsFactors=F,skip=7,nrows=1))
  if(hd%%2 != 0) stop("Data file must have even number of columns!")
  idx <-  matrix(seq_len(hd),ncol=2,byrow=T)
  ret <- list()
  for (i in seq_len(hd/2)) {
    cC1 <- rep("NULL",hd)
    cC2 <- rep("NULL",hd)
    cC1[idx[i,]] <- "character"
    cC2[idx[i,]] <- "numeric"
    tab <- read.table(file,header=T,sep=sep,dec=dec,stringsAsFactors=F,skip=7,colClasses=cC2)
    meta <- read.table(file,header=F,sep=sep,dec=dec,stringsAsFactors=F,skip=0,nrows=5,comment.char="#",colClasses=cC1)
    metadata <- meta[,2]
    names(metadata) <- meta[,1]
    ret[[i]] <- std(tab[,1],tab[,2],desc=as.list(metadata))
  }
  do.call(c.std,ret)
}


#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#:  modify std object
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

tweak.std <- function(x,desc=desc.std(x=x),lmargs=as.list(x[[1]]$lmfit$call[-c(1:3)]))  {
  stopifnot(is.std(x),length(x) == 1,is.list(lmargs),is.list(desc))
  # get data
  a <- x[[1]]$stdata$a
  r <- x[[1]]$stdata$r
  # make new object
  std(a=a,r=r,desc=desc,lmargs=lmargs)
}


#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#:  concatenate
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

c.std <- function(...)  {
  dots <- list(...)
  if(all(sapply(dots,is.std))){
    dots <- lapply(dots,unclass)
    ret <- do.call(c,dots)
    class(ret) <- c("std",class(ret))
    return(ret)
  } else stop("All arguments must have std class.")
}


#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#:  subset
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

'[.std' <- function(x, i, ...) {
structure(NextMethod("["), class = class(x))
}


#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#:  transform from usr to axis coordinates and vise versa
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

y2rry <- function(x) {log(log(100/x))}
rry2y <- function(x) {100/exp(exp(x))}
x2rrx <- function(x) {log(x)}
rrx2x <- function(x) {exp(x)}


#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#:  plot
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

plot.std <-   function(x,
                       type=c("rr","rrdist"),
                       n=30, lgd=T,
                       col=c("#1B9E77", "#D95F02", "#7570B3", "#E7298A",
                             "#66A61E", "#E6AB02", "#A6761D", "#666666"),
                       bx=1,
                       by=5,
                       lgd.x="topleft",
                       ...)
{
  # x ... std oject
  # n ... hint on labels density

  # rounds to nearest b
  mround <- function(x,b=5){
    ord <- floor(log(base=10,x))
    base <- b*10^(ord-1)
    base*round(x/base)
  }

  # pretty labels
  pretty_loglog <- function(x,n,b)  {
    # x ... loglog numbers ~ stdata$Y
    # n ... hint for count
    rge <- range(x,finite=T,na.rm=T)
    unique(mround(rry2y(seq(rge[1],rge[2],length.out=n)),b=b))
  }

  pretty_log <- function(x,n,b)  {
    # x ... log numbers ~ stdata$X
    # n ... hint for count
    rge <- range(x,finite=T,na.rm=T)
    unique(mround(rrx2x(seq(rge[1],rge[2],length.out=n)),b=b))
  }

  stopifnot(is.std(x), is.character(type))

  if (length(x) > length(col)) {
    warning("There is more samples than colors in palette! Forcing rainbow and disabling legend.")
    col <- rainbow(length(x))
    lgd <- F
  }

  dots <- list(...)

  if (is.character(type)) {
    type <- match.arg(type)
  }

  # legend entries
  lgdtext <- sapply(x, FUN=function(x){with(x,paste(Title,Sample,Unit,Date,sep=" | "))})
  bg <- col[seq_len(length(x))]
  # point type
  mypch <- seq_len(length(x))
  lty <- 2
  lty.d <- 3
  lty.u <- 1
  lty.o <- 2
  fill <- NA
  border <- NA

  rr <- expression({
    # arguments
    if(is.null(dots$xlab)) dots$xlab <- "Diameter of mesh aperture [ um ]"
    if(is.null(dots$ylab)) dots$ylab <- "Cumulative Oversize, [ mass % ]"
    if(is.null(dots$main)) dots$main <- "Sieve Test - Rosin - Rammler Plot"
    if(is.null(dots$main)) dots$main <- "Sieve Test - Rosin - Rammler Plot"
    if(is.null(dots$cex.axis)) dots$cex.axis <- 0.7
    if(is.null(dots$las)) dots$las <- 1
    if(is.null(dots$lwd)) dots$lwd <- 1

    # common limits
    limX <- lapply(x, FUN=function(x){x$stdata$X})
    limX <- range(unique(sort(do.call(c,limX))),finite=T,na.rm=T)
    limY <- lapply(x, FUN=function(x){x$stdata$Y})
    limY <- range(unique(sort(do.call(c,limY))),finite=T,na.rm=T)

    # common labels
    # xaxis
    xl <- pretty_log(limX,n,bx)
    xa <- x2rrx(xl)
    # yaxis
    yl <- pretty_loglog(limY,n,by)
    ya <- y2rry(yl)

    # base plot
    arg <- c(list(limX, limY, type='n', axes=FALSE), dots)
    do.call('plot.default',arg)
    box()
    axis(1,at=xa,labels=xl,las=2,cex.axis=dots$cex.axis)
    axis(2,at=ya,labels=yl,las=2,cex.axis=dots$cex.axis)
    abline(col="grey", v=xa, h=ya, lwd=0.5, lty=2)
    abline(col="grey", h=0, lwd=1.5, lty=1)

    for (i in seq_len(length(x))) {
      xi <- x[[i]]

      X <- xi$stdata$X
      Y <- xi$stdata$Y
      arg <- c(list(Y ~ X,bg=col[i],col=col[i],pch=mypch[i]),dots)
      do.call('points',arg)
      usr <- par('usr')
      clip(min(X[is.finite(X)]),max(X),usr[3],usr[4])
      abline(xi$lmfit, col=col[i], lty=lty.o, lwd=dots$lwd)
      do.call("clip", as.list(par('usr')))
    }

  })

  rrdist <- expression({
    # arguments
    if(is.null(dots$xlab)) dots$xlab <- "Diameter of mesh aperture [ um ]"
    if(is.null(dots$ylab)) dots$ylab <- "CDF & CCDF Value (left), PDF Value (right)"
    if(is.null(dots$main)) dots$main <- "Distribution Functions"
    if(is.null(dots$log)) dots$log <- ""
    if(is.null(dots$cex.axis)) dots$cex.axis <- 0.7
    if(is.null(dots$las)) dots$las <- 1
    if(is.null(dots$lwd)) dots$lwd <- 1

    if(n < 1000 | is.null(n)) n <- 1000

    limX <- lapply(x, FUN=function(x){x$stdata$a})
    limX <- range(unique(sort(do.call(c,limX))),finite=T,na.rm=T)
    limY <- c(0,1)
    if (grepl(x=dots$log,pattern="x",ignore.case=T)) limX[1] <- limX[2]/n

    # base plot
    arg <- c(list(limY ~ limX, type='n', axes=FALSE),dots)
    do.call('plot.default',arg)


    i <- which(names(dots) == "log")
    dots <- dots[-i]
    for (i in seq_len(length(x))) {
      xi <- x[[i]]

      # maximalni velikost castice
      xmax <- max(limX)

      # rr parameters
      ex <- xi$RRcoefficients$exponent
      xs <- xi$RRcoefficients$diameter
      yp <- xi$stdata$O/100
      xp <- xi$stdata$a

      # calc
      # d curve
      xx <- seq(0,xmax,length.out=n)
      y_d <- drr(x=xx,ex=ex,xs=xs)
      y_d_pl <- y_d/max(pretty(y_d))
      # oversize, undersize
      y_u <- urr(x=xx,ex=ex,xs=xs)
      y_o <- orr(x=xx,ex=ex,xs=xs)
      
   
      arg <- c(list(y_u ~ xx,type='l',col=col[i],lty=lty.u),dots)
      do.call('lines',arg)
      arg <- c(list(y_o  ~ xx,type='l',col=col[i],lty=lty.o),dots)
      do.call('lines',arg)
      arg <- c(list(y_d_pl  ~ xx,type='l',col=col[i],lty=lty.d),dots)
      do.call('lines',arg)
      arg <- c(list(y=yp,x=xp,col=col[i],pch=mypch[i]),dots)
      do.call('points',arg)
    }
      grid()
      box()
      axis(1,cex.axis=dots$cex.axis)
      axis(2,cex.axis=dots$cex.axis)
      at <- axTicks(side=2)
      lab <- at*max(pretty(y_d))
      axis(4,cex.axis=dots$cex.axis,labels=lab,at=at)

      mypch <- c(mypch,NA,NA,NA)
      lgdtext <- c(lgdtext,"CDF (undersize)","CCDF (oversize vs. measured)","PDF (density)")
      lty <- c(rep(NA,times=length(x)),lty.u,lty.o,lty.d)
  })

  cmd <- get(type)
  eval(cmd)

  if (lgd) {
    legend(x=lgd.x, fill=fill, border=border, pch=mypch, pt.lwd=1, lty=lty, lwd=2 , col=bg, pt.bg=bg, legend=lgdtext,bty='n', cex=1.14*dots$cex.axis)
  }

}


#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#:  summary
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

summary.std <- function(object,...)  {
  su <- lapply(object,FUN=function(x){
    with(x,
{
  data.frame(Date,Title,Sample,Unit,Description,
             RRxs=RRcoefficients$diameter,
             RRex=RRcoefficients$exponent,
             RRms=size$modus,
             RRmn=size$mean,
             RRo90umPPC=size$o90um*100,
             stringsAsFactors=F)
})})
su <- do.call(rbind,su)
rownames(su) <- NULL
return(su)
}


#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#:  test is std?
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

is.std <- function(x)  {
  nms <- c("Description", "Title", "Sample", "Unit", "Date", "stdata", "lmfit", "RRcoefficients","size")
  is_list <- is.list(x)
  all_lists <- all(sapply(x,is.list))
  is_std <- inherits(x, "std")
  all_nms <- all(sapply(x,FUN=function(x){all(names(x) %in% nms)}))
  all_has_lm <- all(sapply(x,FUN=function(x){inherits(x$lmfit, "lm")}))
  all_has_tab <- all(sapply(x,FUN=function(x){inherits(x$stdata, "data.frame")}))
  is_list & all_lists & is_std & all_nms & all_has_lm & all_has_tab
}


#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#:  RR distribution
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

drr <- function(x,ex,xs){(ex/xs)*(x/xs)^(ex-1)*exp(- (x/xs)^ex)}
# drr2 <- function(x,ex,xs){(ex/xs^ex)*(x)^(ex-1)*exp(- (x/xs)^ex)}
urr <- function(x,ex,xs){1-exp(- (x/xs)^ex)}
orr <- function(x,ex,xs){exp(- (x/xs)^ex)}


#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#:  notes
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#   TODO nahodit funkci pro vygenerovani pdf
#   TODO read.std: too many hardcoded arguments
