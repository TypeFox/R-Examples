#' Create a censored time series object of \code{cenTS} class
#'
#' Create a censored time series response object of \code{cenTS} class. Default name of the response is "value",
#'  with the vectors of lower/upper censoring limits denoted by \code{lcl} and \code{ucl} respectively.
#'  The vector of censoring indicators, i.e., \code{ci}, is part of the \code{cenTS} object.
#'  Additional related variables can be stored and provided in the construction function, whose names
#'  are stored in \code{xreg}. All variable values are assumed to be of the same length of and thus
#'  aligned with the censored response time series. \code{cenTS} inherits from \link[xts]{xts}.
#' @param order.by the index vector, must be a vector of time/date.
#' @param value the value vector.
#' @param lcl the vector of lower censoring limits, or a single numeric representing the constant limit.
#'  Default = \code{NULL} indicating no lower limit.
#' @param ucl the vector of upper censoring limits, or a single numeric representing the constant limit.
#'  Default = \code{NULL} indicating no upper limit.
#' @param ci the vector of censoring indicators whose value is -1 (0, 1)
#' if the corresponding response is left censored  (observed, right censored).
#' Default = \code{NULL}, in which case, the function will compute \code{ci} by \code{value}, \code{lcl}
#' and \code{ucl}. If \code{ci} is not \code{NULL}, the function will check the consistency of the data,
#' assuming the observed values less (greater) than or equal to left (right) censoring limits are censored,
#' and are observed otherwise. The function will stop if inconsistent results are found.
#' @param value.name the name of the value, default = "value".
#' @param ... additional variables, must be able to be coerced to a \code{data.frame}.
#' @return a \code{cenTS} object, any censored observation will be replaced by its corresponding censoring limit.
#' @export

#' @examples
#' strDates <- c("2000-01-01", "2000-01-02", "2000-01-03", "2000-01-04", "2000-01-05")
#' ts <- cenTS(value=c(1,-2,1,NA,0),
#'             order.by=as.Date(strDates,"%Y-%m-%d"),
#'             lcl=c(-3,-2,-1,-1,0),
#'             ucl=c(3,2,1,1,1),
#'             x=c(1,1,1,1,1),
#'             y=c(2,2,2,2,2))
#'  print(ts)
#'  print(xreg(ts))
#'  plot(ts)
#'
#' \dontrun{
#' #wrong call, case 1
#'ts <- cenTS(value=c(1,-2,1,NA,0),
#'            order.by=as.Date(strDates,"%Y-%m-%d"),
#'            lcl=c(-3,-2,-1,-1,0),
#'            ucl=c(3,2,1,1,1),
#'            ci =c(-1,-1,1,NA,-1)
#')
#'#wrong call, case 2
#'ts <- cenTS(value=c(1,-2,1,NA,0),
#'            order.by=as.Date(strDates,"%Y-%m-%d"),
#'            lcl=c(-3,-2,-1,-1,0),
#'            ucl=c(3,2,1,1,1),
#'            ci =c(1,-1,1,NA,-1)
#')
#'
#'
#'#wrong call, case 3
#'ts <- cenTS(value=c(1,-2,1,NA,0),
#'            order.by=as.Date(strDates,"%Y-%m-%d"),
#'            lcl=c(-3,-2,-1,-1,0),
#'            ucl=c(3,2,1,1,1),
#'            ci =c(0,-1,0,NA,-1)
#')
#' }
#'
cenTS <- function(value, order.by,
                  lcl = NULL,ucl = NULL,
                  ci = NULL,
                  value.name = "value",
                  ...)
{
  #step0: check ... variables
  xreg <- list(...)
  if(length(xreg)>0)
  {
    xregNames <- names(xreg)

    if(value.name %in% xregNames | "lcl" %in% xregNames | "ucl" %in% xregNames)
      stop(paste("Variable names '",value.name, "', 'lcl', and 'ucl' are reserved for cenTS, but there is at least one of them has(ve) appeared in the list of ... variables."))
    xreg <- data.frame(xreg)
  }
  else
    xreg <- NULL
  #for(x in xreg)
  #{
  #  if(length(x) != length(order.by))
  #    stop("variables in xreg must be of the same length as the time series!")
  #}
  hasCI <- TRUE
  if(is.null(ci))
  {
    ci <- rep(0,length(value))
    hasCI <- FALSE
  }
  ci[!is.finite(value)] <- NA

  if(!is.null(lcl))
  {
    if(length(lcl)==1)
      lcl <- rep(lcl,length(value))
    idx <- is.finite(value) & is.finite(lcl)
    idx2 <- value[idx] <= lcl[idx]
    if(hasCI)
    {
      test <- idx2 != (ci[idx]==-1)
      if(any(test))
        stop("Inconsistency found in data, at index ",seq(1,length(value))[idx][test])
    }
    else
      ci[idx][idx2] <- -1
    value[idx][idx2] <- lcl[idx][idx2]
  }

  if(!is.null(ucl))
  {
    if(length(ucl)==1)
      ucl <- rep(ucl,length(value))
    idx <- is.finite(value) & is.finite(ucl)
    idx2 <- value[idx] >= ucl[idx]
    if(hasCI)
    {
      test <- idx2 != (ci[idx]==1)
      if(any(test))
        stop("Inconsistency found in data, at index ",seq(1,length(value))[idx][test])
    }
    else
      ci[idx][idx2] <- 1
    value[idx][idx2] <- ucl[idx][idx2]
  }

  val <- data.frame(value.name=value)
  names(val) <- c(value.name)

  val$lcl = lcl
  val$ucl = ucl
  val$ci = ci
  if(is.null(xreg))
  {
    ret <- xts::xts(data.frame(val),order.by=order.by)
    attr(ret,"xreg") <- NULL
  }
  else
  {
    ret <- xts::xts(data.frame(val,xreg),order.by=order.by)
    attr(ret,"xreg") <- colnames(xreg)
  }
  attr(ret,"value.name") <- value.name
  attr(ret,"censoring.rate") <- mean(abs(ci[is.finite(ci)]))

  class(ret) <- c('cenTS','xts','zoo')
  invisible(ret)
}

#' Print a \code{cenTS} object
#' @param x a \code{cenTS} object.
#' @param ... not used.
#' @return none.
#' @export
#' @examples
#' strDates <- c("2000-01-01", "2000-01-02", "2000-01-03", "2000-01-04", "2000-01-05")
#' ts <- cenTS(value=c(1,-2,1,NA,0),
#'             order.by=as.Date(strDates,"%Y-%m-%d"),
#'             lcl=c(-3,-2,-1,-1,0),
#'             ucl=c(3,2,1,1,1),
#'             x=c(1,1,1,1,1),
#'             y=c(2,2,2,2,2))
#'  print(ts)
#'

print.cenTS <- function (x,...)
{
  print(xts::as.xts(x))
  cat(paste("\nCensoring rate:",round(attributes(x)$censoring.rate,4),"\n"))
}


#' Return the \code{xreg} part of the \code{cenTS} object
#' @param object a \code{cenTS} object.
#' @return the list in \code{xreg}.
#' @seealso \code{\link{cenTS}}.
#' @export
#' @examples
#' strDates <- c("2000-01-01", "2000-01-02", "2000-01-03", "2000-01-04", "2000-01-05")
#' ts <- cenTS(value=c(1,-2,1,NA,0),
#'             order.by=as.Date(strDates,"%Y-%m-%d"),
#'             lcl=c(-3,-2,-1,-1,0),
#'             ucl=c(3,2,1,1,1),
#'             x=c(1,1,1,1,1),
#'             y=c(2,2,2,2,2))
#'  xreg(ts)
#'

xreg <- function(object) UseMethod("xreg")

#' Return the \code{xreg} part of the \code{cenTS} object
#' @param object a \code{cenTS} object.
#' @return the list in \code{xreg}.
#' @seealso \code{\link{cenTS}}.
#' @export
#' @examples
#' strDates <- c("2000-01-01", "2000-01-02", "2000-01-03", "2000-01-04", "2000-01-05")
#' ts <- cenTS(value=c(1,-2,1,NA,0),
#'             order.by=as.Date(strDates,"%Y-%m-%d"),
#'             lcl=c(-3,-2,-1,-1,0),
#'             ucl=c(3,2,1,1,1),
#'             x=c(1,1,1,1,1),
#'             y=c(2,2,2,2,2))
#'  xreg(ts)
xreg.cenTS <- function(object)
{
  if(!is.null(attributes(object)$xreg))
    xts::as.xts(object[,attributes(object)$xreg])
  else
    NULL
}


#' Plot a \code{cenTS} object
#' @param x a \code{cenTS} object.
#' @param type,auto.grid,major.ticks,minor.ticks,major.format,bar.col,candle.col,ann,axes,ylim,main,...
#' standard parameters to control the plot.
#' @seealso \code{\link[xts]{plot.xts}}.
#' @export
#' @examples
#' strDates <- c("2000-01-01", "2000-01-02", "2000-01-03", "2000-01-04", "2000-01-05")
#' ts <- cenTS(value=c(1,-2,1,NA,0),
#'             order.by=as.Date(strDates,"%Y-%m-%d"),
#'             lcl=c(-3,-2,-1,-1,0),
#'             ucl=c(3,2,1,1,1),
#'             x=c(1,1,1,1,1),
#'             y=c(2,2,2,2,2))
#'  plot(ts)
plot.cenTS <- function(x, type = "l", auto.grid = TRUE, major.ticks = "auto",
    minor.ticks = TRUE, major.format = TRUE, bar.col = "grey",
    candle.col = "white", ann = TRUE, axes = TRUE,ylim=NULL,main=NULL, ...)
{

    value.name <- attributes(x)$value.name
    series.title <- series.title <- deparse(substitute(x))
    if(value.name != "value")
      series.title <- value.name
    if(!is.null(main))
      series.title <- main

    ep <- xts::axTicksByTime(x, major.ticks, format.labels = major.format)
    otype <- type

    xycoords <- grDevices::xy.coords(xts::.index(x), x[,value.name])
    ylim0 <- range(zoo::coredata(x)[,intersect(colnames(zoo::coredata(x)),c(value.name,"lcl","ucl"))],na.rm=TRUE,finite=TRUE)
    if(is.null(ylim))
      ylim <- ylim0
    else
    {
      ylim[1] <- min(ylim0[1],ylim[1])
      ylim[2] <- max(ylim0[2],ylim[2])
    }

    graphics::plot(xycoords$x, xycoords$y, type = type, axes = FALSE, ann = FALSE, ylim=ylim,...)


    if(!is.null(x$lcl) && any(is.finite(x$lcl)))
    {
      graphics::lines(xycoords$x,x$lcl,col="red",lty=3)
      idx <- is.finite(x[,value.name]) & is.finite(x$lcl)
      idx2 <- x[,value.name][idx] <= x$lcl[idx]
      graphics::points(xycoords$x[idx][idx2],x$lcl[idx][idx2],pch=2)
    }

    if(!is.null(x$ucl) && any(is.finite(x$ucl)))
    {
      graphics::lines(xycoords$x,x$ucl,col="red",lty=5)
      idx <- is.finite(x[,value.name]) & is.finite(x$ucl)
      idx2 <- x[,value.name][idx] >= x$ucl[idx]
      graphics::points(xycoords$x[idx][idx2],x$ucl[idx][idx2],pch=6)
    }

    if (auto.grid) {
        graphics::abline(v = xycoords$x[ep], col = "grey", lty = 4)
        graphics::grid(NA, NULL)
    }
    dots <- list(...)
    if (axes) {
        if (minor.ticks)
            graphics::axis(1, at = xycoords$x, labels = FALSE, col = "#BBBBBB",
                ...)
        graphics::axis(1, at = xycoords$x[ep], labels = names(ep), las = 1,
            lwd = 1, mgp = c(3, 2, 0), ...)
        graphics::axis(2, ...)
    }
    graphics::box()
    if (!"main" %in% names(dots))
        graphics::title(main = series.title)
    do.call("title", list(...))
}
