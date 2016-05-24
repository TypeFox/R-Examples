freq2d <- function(x, ...)
{
  UseMethod("freq2d")
}


freq2d.formula <- function(formula, data, subset, ...)
{
  m <- match.call(expand.dots=FALSE)
  if(is.matrix(eval(m$data,parent.frame())))
    m$data <- as.data.frame(data)
  m$... <- NULL
  m[[1L]] <- as.name("model.frame")
  mf <- eval(m, parent.frame())

  freq2d.default(mf[2:1], ...)
}


freq2d.default <- function(x, y, n=20, pad=0, layout=1, print=TRUE, dnn=NULL, ...)
{
  method <- match.arg(as.character(layout), c("1","2","3"))

  dnn <- if(!is.null(dnn)) rep(dnn,length.out=2) else NULL
  xname <- dnn[1]
  yname <- dnn[2]

  ## 1  Extract data
  if(is.matrix(x))
    x <- as.data.frame(x)
  if(is.list(x))  # data.frame or list
  {
    xname <- if(is.null(xname)) names(x)[1] else xname
    yname <- if(is.null(yname)) names(x)[2] else yname
    y <- x[[2]]
    x <- x[[1]]
  }

  ## 2  Create grid
  n <- rep(n, length.out=2)
  xmid <- pretty(x, n=n[1])
  xstep <- diff(xmid)[1]
  xgrid <- c(xmid-0.5*xstep, max(xmid)+0.5*xstep)
  ymid <- pretty(y, n=n[2])
  ystep <- diff(ymid)[1]
  ygrid <- c(ymid-0.5*ystep, max(ymid)+0.5*ystep)

  ## 3  Map data on grid
  xfac <- cut(x, xgrid, include.lowest=TRUE, labels=format(xmid))
  if(is.null(xname))
    xname <- deparse(substitute(x))
  yfac <- cut(y, ygrid, include.lowest=TRUE, labels=format(ymid))
  if(is.null(yname))
    yname <- deparse(substitute(y))
  z <- table(xfac, yfac, dnn=c(xname,yname))

  ## 4  Remove existing edges with only zeros
  z <- z[cumsum(rowSums(z))>0, cumsum(colSums(z))>0]
  z <- z[rev(cumsum(rev(rowSums(z))))>0, rev(cumsum(rev(colSums(z))))>0]

  ## 5  Add edges with only zeros
  for(i in seq_len(pad))
  {
    tmp <- cbind(0, rbind(0, z, 0), 0)
    rownames(tmp)[c(1,nrow(tmp))] <- as.numeric(rownames(z)[c(1,nrow(z))]) + c(-xstep,xstep)
    colnames(tmp)[c(1,ncol(tmp))] <- as.numeric(colnames(z)[c(1,ncol(z))]) + c(-xstep,xstep)
    names(dimnames(tmp)) <- names(dimnames(z))
    z <- tmp
  }

  ## 5  Prepare output
  xnum <- as.numeric(rownames(z))
  ynum <- as.numeric(colnames(z))
  if(layout == 1)
  {
    output <- t(z)[ncol(z):1,]
    if(print)
    {
      print.table(output, zero.print=".")
      return(invisible(output))
    }
    else
    {
      return(output)
    }
  }
  else if(layout == 2)
  {
    output <- list(x=xnum, y=ynum, z=z)
    return(output)
  }
  else  # layout 3
  {
    output <- data.frame(x=rep(xnum,length(ynum)), y=rep(ynum,each=length(xnum)), z=c(z))
    names(output) <- make.names(c(xname,yname,"Freq"), unique=TRUE)
    return(output)
  }
}
