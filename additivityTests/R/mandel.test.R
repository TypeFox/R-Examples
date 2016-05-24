#' Mandel Additivity Test
#' 
#' Test for an interaction in two-way ANOVA table by the Mandel test.
#'
#' @param data data matrix
#' @param alpha level of the test
#' @param critical.value result of \code{\link{critical.values}} function, see \code{Details}
#'
#' @return A list with class "\code{aTest}" containing the following components: 
#' test statistics \code{stat}, critical value \code{critical.value} and the result of 
#' the test \code{result}, i.e. whether the additivity hypothesis has been rejected.
#'
#' @details The critical value can be computed in advance and given in the parameter \code{critical value}. 
#' If not a function  \code{\link{critical.values}} is called to do that.
#' 
#' @references Mandel, J.: Non-additivity in Two-way Analysis of Variance,
#' \emph{Journal of the American Statistical Association} \bold{56},
#' pp. 878--888, 1961.
#' 
#' @seealso \code{\link{tukey.test}}, \code{\link{mtukey.test}}, \code{\link{johnson.graybill.test}}, 
#' \code{\link{lbi.test}}, \code{\link{tusell.test}}
#' 
#' @keywords htest
#'
#' @export
#' 
#' @examples
#' data(Boik)
#' mandel.test(Boik)

`mandel.test` <-
function(data,alpha=0.05,critical.value=NA) {
# Mandel test of additivity in two-way ANOVA
# In rows is factor A (fixed factor), in columns factor B (random)

  a=nrow(data) # number of levels of factor A
  b=ncol(data)	# number of levels of fator B 
  d.f1=a-1 # df for SS.mandel
  d.f2=(a-1)*(b-2) #df for SS.resid
  if (is.na(critical.value)) critical.value=qf(1-alpha,,df1=d.f1,df2=d.f2)

  yMEAN=mean(data) # grand mean
  A.hat=apply(data,1,mean)-yMEAN # deviations of the row means from the grand mean
  B.hat=apply(data,2,mean)-yMEAN # deviations of the column means from the grand mean
  ss.col=a*sum(B.hat^2) # SS of columns
  ss.row=b*sum(A.hat^2) # SS of rows
  kappa=(data%*%B.hat)/(ss.col/a)  #(sum_j^B y_{ij}*b.hat_j)/(sum_j b.hat_j^2)
  ss.total=sum(apply(data^2,1,sum))  
       #SS total, in accordance with Mandel's (1961) ANOVA layout; (ss.total-ss.mean) is the usual way of expressing it as (sum_i sum_j (y_{ij}-y..2)
  ss.mean=a*b*yMEAN^2 #SS for mean, in accordance with Mandel's (1961) ANOVA layout 
  ss.mandel=sum((kappa-1)^2)*(ss.col/a) #SS for Interaction that is due to differences in the Slopes 
  ss.resid=ss.total-ss.mean-ss.row-ss.col-ss.mandel # Residual SS
  test.stat=(ss.mandel/d.f1)/(ss.resid/d.f2) # The F statistic as MS.mandel/MS.resid

out<-list(result=test.stat>critical.value,stat=test.stat,critical.value=critical.value,alpha=alpha,name="Mandel test") # nezamitame aditivitu

  class(out)<-"aTest"
  return(out)
}

