col2hex <- function(cname)
  {
    colMat <- col2rgb(cname)
    rgb(
        red=colMat[1,]/255,
        green=colMat[2,]/255,
        blue=colMat[3,]/255
        )
  }
