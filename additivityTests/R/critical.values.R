#' Critical Values for the Johnson-Graybill, LBI and Tusell tests
#' 
#' Compute the critical values by performing N simulation.
#'
#' @param a number of rows
#' @param b number of columns
#' @param N number of simulations
#' @param alpha level(s) of the test
#'
#' @return A list containing three components: critical values for Johnson-Graybill, LBI and Tusell tests, respectively.
#'
#' @seealso \code{\link{johnson.graybill.test}}, \code{\link{lbi.test}}, \code{\link{tusell.test}}
#' 
#' @keywords htest
#'
#' @export
#' 
#' @examples
#' data(Boik)
#' critical.values(nrow(Boik), ncol(Boik), 0.01) 

`critical.values` <-
function(a, b, N=100000, alpha=0.05)  
{

  t1<-t2<-t3<-NULL
  p<-a-1
  q<-b-1

  for (i in 1:N)
  {

    Y<-matrix(rnorm(a*b),a,b)
    R<-Y-rep(apply(Y,1,mean),b)-rep(apply(Y,2,mean),each=a)+rep(mean(Y),a*b)
    S<-R %*% t(R)
    vl.cisla<-eigen(S / sum(diag(S)),only.values = TRUE)$values

    # Johnson & Graybill test
    t1<-c(t1,vl.cisla[1])
    # LBI test 
    t2<-c(t2,sum(vl.cisla^2))
    # Tusell test
    t3<-c(t3,prod(vl.cisla[1:p]))
  }

  return(list(t1=quantile(t1,1-alpha),t2=quantile(t2,1-alpha),t3=quantile(t3,alpha),alpha=alpha))
}
