#' @title Compute standardized tracheidograms
#' @description
#' This function computes standardized tracheidograms from raw tracheidograms.
#' @param traq  a raw tracheidogram (or a list of raw tracheidograms)
#' @param method a \code{string} defining the \code{method} to be used. Possible values are c("nCells","kCells","relPos"). The default method is "nCells".
#' @param k a integer to determine the number of cells of the standardized tracheidogram when \code{method} is "kCells"
#' @details
#' This function takes as input raw tracheidograms (obtained using the \code{getTrac} function) and standardizes them using 3 different methods. The first method ("nCells") standardizes rays of a given ring using the mean number of cells, allowing that different rings have different number of cells. The method "kCells" normalizes tracheidograms to a constant number of tracheids (\code{k}). The "relPos" method standardizes the tracheidogram based on the relative position of each tracheid inside the tree ring.
#' @return the function \code{tracheider} returns a list with the following elements:
#' @return \code{LD} ordered sequence of lumen diameters.
#' @return \code{CWT} ordered sequence of radial cell wall thikness.
#' @return \code{LWratio} ordered sequence of LD/CWT ratio.
#' @export
#' @importFrom tgram tgram
#' @examples
#'
#' ## Not run:
#' y2010 <- getTrac(tch$'y2010', scale=0.169)
#' y2013 <- getTrac(tch$'y2013', scale=0.169)
#'
#' ## nCells
#' y2010n <- tracheider(y2010, method = "nCells")
#' y2013n <- tracheider(y2013, method = "nCells")
#' plot(rowMeans(as.data.frame(y2010n$LD)), type="l", xlim=c(0,86),
#'      ylim=c(0,60), yaxs="i", xaxs="i", xlab="Number of tracheid",
#'      ylab=expression(paste("LD (", mu,"m)")), col=2, lwd=2)
#' lines(rowMeans(as.data.frame(y2013n$LD)), col="blue", lwd=2)
#' legend("topright",lty=1,lwd=2,col = c(2,4), legend=c("2010 ","2013 "),
#'         text.col = c(2,4), box.col = "#00000000", bg="#00000000")
#'
#' ## k = 53
#' TCH <- lapply(tch, getTrac, scale=0.169)
#' TCHn53 <- lapply(TCH, tracheider,method = "kCells", k=53)
#'
#' plot(rowMeans(as.data.frame(TCHn53$'y2010'$LD)), type="l",
#'      xlab="Number of tracheid", ylab=expression(paste("LD (", mu,"m)")),
#'      xlim=c(0,54), ylim=c(0,60), col=2, lwd=2, yaxs="i", xaxs="i")
#' lines(rowMeans(as.data.frame(TCHn53$'y2013'$LD)), col="blue", lwd=2)
#' legend("topright",lty=1,lwd=2,col = c(2,4),legend=c("2010 ","2013 "),
#'      text.col = c(2,4), box.col = "#00000000", bg="#00000000")
#'
#' ## Relative position
#' TCH <- lapply(tch, getTrac, scale=0.169)
#' TCHrelPos <- lapply(TCH, tracheider, method ="relPos")
#' plot(rowMeans(as.data.frame(TCHrelPos$'y2010'$LD)),
#'      type="l", xlim=c(0,101), ylim=c(0,60), col=2, lwd=2,
#'      xlab="Number of tracheid", ylab=expression(paste("LD (", mu,"m)")),
#'      yaxs="i", xaxs="i")
#' lines(rowMeans(as.data.frame(TCHrelPos$'y2013'$LD)),
#'       col="blue", lwd=2)
#' legend("topright",lty=1,lwd=2,col = c(2,4), legend=c("2010 ","2013 "),
#'        text.col = c(2,4), box.col = "#00000000", bg="#00000000")
#' ## End(not run)


tracheider = function (traq, method = c("nCells","kCells","relPos"), k = 20) {
  method <- match.arg(method)

  getN = function (x) {
    n = length(x)
    out = rep(NA,n)
    for (i in 1:n)
      out[i] <- x[[i]]$n
    round(mean(out))
  }

#   #stdN versus standz
#   getTrac(tch$y2010$ray1, scale = 0.169) -> y2010r1
#   y2010r1$LD -> y
#   n = 40
#   plot(1:13,y, ann = FALSE, axes = FALSE)
#   lines(seq(1,length(y),(length(y) - 1) / (n - 1)),standz(y,G = n),
#         col = 2)
#   lines(seq(1,length(y),(length(y) - 1) / (n - 1)),
#         approx(1:length(y),y,seq(1,length(y),(length(
#           y) - 1) / (n - 1)))$y,
#         col = 4)
#   axis(1, at = 1:13)
#   axis(2)


  stdN = function (y, n)
    approx(1:length(y),y,seq(1,length(y),(length(y) - 1) / (n - 1)))$y

  STD = function (x, n = NULL) {
    if (is.null(n))
      n = getN(x)
    LDstd <- lapply(x, function (x)
      stdN(x$LD, n = n))
    CWTstd <- lapply(x, function (x)
      stdN(x$CWT, n = n))
    LWratio.std <- lapply(x, function (x)
      stdN(x$LWratio, n = n))
    return(list(
      "LD" = LDstd,"CWT" = CWTstd, "LWratio" = LWratio.std
    ))
  }

  relPos =  function (x) {
    stdPos =  function (x, y) {
      c(0,x$pos$x / x$pos$RingWidth, 1) * 100 -> x
      linearInterp(x, c(y[1],y,y[length(y)]))$y
    }
    LD = lapply(x, function (x)
      stdPos(x,x$LD))
    CWT = lapply(x, function (x)
      stdPos(x,x$CWT))
    LWratio = lapply(x, function (x)
      stdPos(x,x$LWratio))
    return(list(
      "LD" = LD, "CWT" = CWT, "LWratio" = LWratio
    ))

  }

  switch(
    method,
    "nCells" = return(STD(traq)),
    "kCells" = return(STD(traq, n = k)),
    "relPos" = return(relPos(traq))
  )

}
