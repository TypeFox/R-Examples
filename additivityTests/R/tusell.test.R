#' Tusell Additivity Test
#' 
#' Test for an interaction in two-way ANOVA table by the Tusell test.
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
#' @references Tusell, F.: Testing for Interaction in Two-way ANOVA Tables with no Replication, 
#' \emph{Computational Statistics \& Data Analysis} \bold{10}, pp. 29--45, 1990
#' 
#' @seealso \code{\link{tukey.test}}, \code{\link{mtukey.test}}, \code{\link{mandel.test}}, 
#' \code{\link{lbi.test}}, \code{\link{johnson.graybill.test}}
#' 
#' @keywords htest
#'
#' @export
#' 
#' @examples
#' data(Boik)
#' tusell.test(Boik)

`tusell.test` <-
function(Y,alpha=0.05,critical.value=NA,Nsim=1000)
{

  if (nrow(Y)>ncol(Y)) Y<-t(Y)
  if (is.na(critical.value)) critical.value<-critical.values(nrow(Y),ncol(Y),Nsim,alpha)$t3

  a<-nrow(Y)
  b<-ncol(Y)
  p<-a-1
  q<-b-1

  R<-Y-rep(apply(Y,1,mean),b)-rep(apply(Y,2,mean),each=a)+rep(mean(Y),a*b)
  S<-R %*% t(R)
  vl.cisla<-eigen(S / sum(diag(S)),only.values = TRUE)$values

  if (prod(vl.cisla[1:p])<critical.value)  out<-list(result=TRUE,stat=prod(vl.cisla[1:p]),critical.value=critical.value,alpha=alpha,name="Tusell test") # zamitame aditivitu
                             else  out<-list(result=FALSE,stat=prod(vl.cisla[1:p]),critical.value=critical.value,alpha=alpha,name="Tusell test") # nezamitame aditivitu

  class(out)<-"aTest"
  return(out)

}

