#' Johnson and Graybill Additivity Test
#' 
#' Test for an interaction in two-way ANOVA table by the Johnson-Graybill test.
#'
#' @param Y data matrix
#' @param alpha level of the test
#' @param critical.value result of \code{\link{critical.values}} function, see \code{Details}
#' @param Nsim number of simulations to be used for a critical value estimation
#'
#' @return A list with class "\code{aTest}" containing the following components: 
#' test statistics \code{stat}, critical value \code{critical.value} and the result of 
#' the test \code{result}, i.e. whether the additivity hypothesis has been rejected.
#'
#' @details The critical value can be computed in advance and given in the parameter \code{critical value}. 
#' If not a function  \code{\link{critical.values}} is called to do that.
#' 
#' @references Johnson, D.E. and Graybill, F.A.: An analysis of a two-way model with interaction and no replication, 
#' \emph{Journal of the American Statistical Association} \bold{67}, pp.  862--868, 1972.
#' 
#' @seealso \code{\link{tukey.test}}, \code{\link{mtukey.test}}, \code{\link{mandel.test}}, 
#' \code{\link{lbi.test}}, \code{\link{tusell.test}}
#' 
#' @keywords htest
#'
#' @export
#' 
#' @examples
#' data(Boik)
#' johnson.graybill.test(Boik)

`johnson.graybill.test` <-
function(Y, alpha=0.05, critical.value=NA, Nsim=1000)
{

  if (nrow(Y)>ncol(Y)) Y<-t(Y)
  if (is.na(critical.value)) critical.value<-critical.values(nrow(Y),ncol(Y),Nsim,alpha)$t1

  a<-nrow(Y)
  b<-ncol(Y)
  p<-a-1
  q<-b-1

  R<-Y-rep(apply(Y,1,mean),b)-rep(apply(Y,2,mean),each=a)+rep(mean(Y),a*b)
  S<-R %*% t(R)
  vl.cisla<-eigen(S / sum(diag(S)),only.values = TRUE)$values

  if (vl.cisla[1]>critical.value) out<-list(result=TRUE,stat=vl.cisla[1],critical.value=critical.value,alpha=alpha,name="Johnson and Graybill test") # zamitame aditivitu
                             else out<-list(result=FALSE,stat=vl.cisla[1],critical.value=critical.value,alpha=alpha,name="Johnson and Graybill test")

  class(out)<-"aTest"
  return(out)

}

