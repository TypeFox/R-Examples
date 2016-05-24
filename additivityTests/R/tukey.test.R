#' Tukey Additivity Test
#' 
#' Test for an interaction in two-way ANOVA table by the Tukey test.
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
#' @references Tukey, J.W.: One Degree of Freedom for Non-additivity, \emph{Biometrics} \bold{5}, 
#' pp. 232--242, 1949.
#' 
#' @seealso \code{\link{tusell.test}}, \code{\link{mtukey.test}}, \code{\link{mandel.test}}, 
#' \code{\link{lbi.test}}, \code{\link{johnson.graybill.test}}
#' 
#' @keywords htest
#'
#' @export
#' 
#' @examples
#' data(Boik)
#' tukey.test(Boik)

`tukey.test` <-
function(data, alpha=0.05, critical.value=NA) {
# Mandel test of additivity in two-way ANOVA
# In rows is factor A (fixed factor), in columns factor B (random)

  a<-nrow(data) # number of levels of factor A
  b<-ncol(data) # number of levels of factor B	

  d.f1=1 # df for SS.mandel
  d.f2=(a-1)*(b-1)-1 #df for SS.resid
  if (is.na(critical.value)) critical.value=qf(1-alpha,df1=d.f1,df2=d.f2)

  yMEAN<-mean(data) # grand mean
  A.hat<-apply(data,1,mean)-yMEAN # deviations of the row means from the grand mean
  B.hat<-apply(data,2,mean)-yMEAN # deviations of the row means from the grand mean
  ss.col<-a*sum(B.hat^2) # SS of columns	
  ss.row<-b*sum(A.hat^2) # SS of rows
  ss.tukey<-((sum(A.hat*(data%*%B.hat)/(ss.col/a)))^2)/(ss.row/b)*(ss.col/a) 

  test.stat=ss.tukey*((a-1)*(b-1)-1)/(sum(apply(data^2,1,sum))-a*b*yMEAN^2-ss.row-ss.col-ss.tukey) # The F statistic
 out<-list(result=test.stat>critical.value,stat=test.stat,critical.value=critical.value,alpha=alpha,name="Tukey test") # nezamitame aditivitu

  class(out)<-"aTest"
  return(out)
}

