#' @importFrom stats approx
linearInterp = function (x, y, start=1, end=100){
  approx(x, y, xout=start:end)
}
