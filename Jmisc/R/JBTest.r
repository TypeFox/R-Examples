#' Return the p Value of Jarque Bera test. The Jarque Bera test  test the null hypothesis that the data are from a normal distribution. 
#' @name JBTest
#' @aliases JBTest
#' @title p Value of Jarque Bera test
#' @param x data
#' @return p Value of Jarque Bera test
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
#' @examples        
#' JBTest(rnorm(50))
#' JBTest(rt(50,3))
#' 
#' n=100
#' # size
#' mean(replicate(n,JBTest(rnorm(100)))<0.05)
#' 
#' # power
#' mean(replicate(n,JBTest(rt(100,3)))<0.05)
JBTest <-
function(x){
    m1<-mean(x)
    m2<-mean((x-m1)^2)
    m3<-mean((x-m1)^3)
    m4<-mean((x-m1)^4)
    s<-m3/m2^(3/2)
    k<-m4/m2^2-3
    jb_test=length(x)/6*(s^2+k^2/4)
    1 - pchisq(jb_test, df = 2)
}

