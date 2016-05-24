if(getRversion() >= '2.15.1') utils::globalVariables(c("item"))

#' Functions to find out constants of the relative range distribution.
#' 
#' These functions compute the constants d2, d3 and c4 to get estimators of the
#' standard deviation to set control limits.
#' 
#' @name ss.cc.constants
#' @aliases ss.cc.getd2 ss.cc.getd3 ss.cc.getc4
#' @usage ss.cc.getd2(n = NA)
#' @usage ss.cc.getd3(n = NA)
#' @usage ss.cc.getc4(n = NA)
#' @param n Sample size
#' @return A numeric value for the constant.
#' 
#' @references 
#' Cano, Emilio L., Moguerza, Javier M. and Redchuk, Andres. 2012.
#' \emph{Six Sigma with {R}. Statistical Engineering for Process
#'   Improvement}, Use R!, vol. 36. Springer, New York.
#'   \url{http://www.springer.com/statistics/book/978-1-4614-3651-5}.
#' 
#' @seealso ss.cc
#' 
#' @author EL Cano
#' @examples 
#' ss.cc.getd2(20)
#' ss.cc.getd3(20)
#' ss.cc.getc4(20)
#' @keywords control charts constants
#' @export ss.cc.getd2 ss.cc.getd3 ss.cc.getc4
ss.cc.getc4 <- function(n = NA){
  if (is.na(n) | n < 2 | abs(n-round(n))!=0){
    stop("Invalid sample size")
  }
  c4=sqrt(2/(n-1)) * ((gamma(n/2))/(gamma((n-1)/2)))
  names(c4)=c("c4")
  return(c4)
}

ss.cc.getd2 <- function (n = NA){
  if (is.na(n) | n < 2 | abs(n-round(n))!=0){
    stop("Invalid sample size")
  }
  f <- function(x){
    1 - ptukey(x, n, Inf)
  }
  d2 <- integrate(f, 0, Inf)
  if (d2$abs.error > 0.001) 
    warning("The absolute error after numerical integration 
            is greater than 0.001")
  d2 <- d2$value
  names(d2) <- c("d2")
  return(d2)
}

ss.cc.getd3 <- function (n = NA){
  if (is.na(n) | n < 2 | abs(n-round(n))!=0){
    stop("Invalid sample size")
  }
  f <- function (x){
    x * (1 - ptukey(x, n, Inf))
  }
  d3 <- integrate(f, 0, Inf)
  if (d3$abs.error > 0.001) 
    warning("The absolute error after numerical integration
            is greater than 0.001")
  d3 <- 2 * d3$value
  this.d2 <- ss.cc.getd2(n)
  d3 <- sqrt(d3 - this.d2^2)
  names(d3) <- c("d3")
  return(d3)
}


#' Control Charts
#' 
#' Plot control charts
#' 
#' If control limits are provided, \code{cdata} is dismissed and a message is 
#' shown. If there are no control limits nor controlled data, the limits are 
#' computed using \code{data}.
#' \cr
#' Supported types of control charts:
#' \itemize{
#' \item{mr}{Moving Range}
#' }
#' @note 
#' We have created this function since the \code{qAnalyst} package has been
#' removed from \code{CRAN}, and it was used in the "Six Sigma with R" book to
#' plot moving average control charts.
#' 
#' 
#' @param type  Type of chart (see details)
#' @param data data.frame with the process data.
#' @param cdata Vector with the controlled process data to compute limits.
#' @param CTQ Name of the column in the data.frame containing the CTQ.
#' @param groups  Name of the column in the data.frame containing the groups.
#' @param climits Limits of the controlled process. It should contain three 
#' ordered values: lower limit, center line and upper limit. 
#' @param nsigmas Number of sigmas to compute the limits from the center line
#' (default is 3)
#' 
#' @return 
#' A plot with the control chart, and a list with the following elements:
#' \item{LCL}{Lower Control Limit}
#' \item{CL}{Center Line}
#' \item{UCL}{Upper Control Limit}
#' \item{phase}{II when cdata or climits are provided. I otherwise.}
#' \item{out}{Out of control points}
#' 
#' @references 
#' Cano, Emilio L., Moguerza, Javier M. and Redchuk, Andres. 2012.
#' \emph{Six Sigma with {R}. Statistical Engineering for Process
#'   Improvement}, Use R!, vol. 36. Springer, New York.
#'   \url{http://www.springer.com/statistics/book/978-1-4614-3651-5}.
#' 
#' @seealso \code{\link{ss.cc.constants}}
#' @author EL Cano
#' @examples 
#' ss.cc("mr", ss.data.pb1, CTQ = "pb.humidity")
#' testout <- ss.data.pb1
#' testout[31,] <- list(31,17)
#' ss.cc("mr", testout, CTQ = "pb.humidity")
#' 
#' @export
ss.cc <- function(type, data, cdata, CTQ = names(data)[1], groups,
    climits, nsigmas = 3){
  # Input control
  if (!is.data.frame(data)){
    stop("Please provide a data frame as input data")
  }
  if (!is.numeric(as.data.frame(data)[, eval(CTQ)])){
    stop("Incorrect data.frame column")
  }
  if (length(data[, eval(CTQ)]) < 2){
    stop("Not enough data to plot")
  }
  
  if (!missing(climits) && (length(climits) != 3 || (1:length(climits) != order(climits)))){
    stop("Incorrect control limits.")
  }
  if (!missing(cdata) && !missing(climits)){
    message("The control limits in argument 'climits' are used. 'cdata'
            argument is ignored.")
  }
  
  # If control limits are provided, they are directly used
  # otherwise the appropriate data source is selected to compute them
  if (!missing(climits)){
    LCL <- climits[1]
    CL <- climits[2] 
    UCL <- climits[3]
    phase <- "II"
  } else if (!missing(cdata)){
    data.lim <- cdata
    phase <- "II"
  } else {
    data.lim <- data
    phase <- "I"
  }
  # Select type of control chart. Stop if not supported. 
  if (type == "mr"){
    # Compute lines (center and limits)
    if (missing(climits)){
      datalim.vector <- data.lim[, eval(CTQ)]
      MRlim <- sapply(1:(length(datalim.vector) - 1), function(x){
          abs(datalim.vector[x + 1] - datalim.vector[x])
        })
#      avgMR <- mean(MR)
      CL <- mean(MRlim)
      LCL <- CL * (1 - nsigmas * (ss.cc.getd3(2) / ss.cc.getd2(2)))
      UCL <- CL * (1 + nsigmas * (ss.cc.getd3(2) / ss.cc.getd2(2)))
    }
    if (LCL < 0) LCL <- 0 # Always matters if nsigmas > 1.32 (usual case)
    # Get points to plot
    data.vector <- data[, eval(CTQ)]
    MR <- sapply(1:(length(data.vector) - 1), function(x){
          abs(data.vector[x + 1] - data.vector[x])
        })
    gdata <- data.frame(item = 1:length(MR), 
        MR, out = MR > UCL | MR < LCL)
    outData <- subset(gdata, out == TRUE)
    # Plot chart using ggplot2 library
    ccPlot <- ggplot(data = gdata, 
            aes(x = item, y = MR)) + 
        geom_point() +
        geom_line() +
        geom_hline(yintercept = c(LCL, UCL), size = 1) +
        geom_hline(yintercept = CL, linetype = 2) +
        ggtitle("Moving Range Control Chart") +
        scale_y_continuous(name = "Moving Range")
    if (nrow(outData) > 0){
      ccPlot <- ccPlot + geom_point(x = outData$item, 
          y = outData$MR, 
          col = "red",
          size = 2.5)
    }
    
  } else{
    stop("Unsupported type of control chart.")
  }
  
  out <- list(
      LCL = as.vector(LCL), 
      CL = as.vector(CL), 
      UCL = as.vector(UCL), 
      phase = phase,
      out = as.vector(outData$item))
  print(ccPlot)
  cat(paste("Phase", out[["phase"]], "limits:\n"))
  print(unlist(out[1:3]))
  cat("\nOut of control Moving Range:\n")
  if (length(out[["out"]]) > 0){
    print(out[["out"]])
  } else {
    cat("None\n")
  }
  invisible(out)
}




