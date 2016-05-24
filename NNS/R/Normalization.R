#' VN Normalization
#'
#' Normalizes a matrix of variables based on nonlinear scaling normalization method.
#' @param A Matrix of variables.
#' @param chart.type  Defaults to NULL.  'l' for line, 'b' for boxplot
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{http://amzn.com/1490523995}
#' @examples
#' set.seed(123)
#' x<-rnorm(100); y<-rnorm(100)
#' A<-cbind(x,y)
#' \dontrun{VN.norm(A, chart.type='l')}
#' @export

VN.norm <- function(A,chart.type=NULL) {
  m  <- colMeans(A)
  RG <- m %o% (1/m)
  scales <- colMeans(RG * abs(cor(A)))
  A_Normalized <- t(t(A) * scales)

  n <- ncol(A)
  i <- seq_len(n)
  labels <- c(colnames(A),
              paste0(colnames(A)," Normalized"))


if(!is.null(chart.type)){
  if(chart.type== 'b' ){
  par(mar = c(10,4,3,2) + 0.1)
  boxplot(cbind(A, A_Normalized),
          las = 2, names = labels,
          col = c(rep("grey", n), rainbow(n)))

  }

  if(chart.type== 'l' ){

    par(mar = c(3,3,3,2))

  par(mfrow=c(2,1))

  matplot(A,type = 'l',col=rainbow(n),ylab='')
   legend('topleft', inset=c(0,0),c(colnames(A)),lty=1,col=c(rainbow(n)))

  matplot(A_Normalized,type = 'l',col=rainbow(n),ylab='')
  legend('topleft',c(paste0(colnames(A)," Normalized")),lty=1,col=c(rainbow(n)))
  }}

  par(mfrow=c(1,1))
  return(A_Normalized)

}
