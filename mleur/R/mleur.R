mleur <-
function(y, type=c("p", "n")){
#implements mle unit root test using the results for the response-surface
#SHARCNET using 221 CPUs to compute.
#The settings are M=200, N=200,000
#sample sizes chosen as
# n=c(seq(20,100,5),seq(120,300,20),seq(350,500, 50), seq(600,1000,100))
  type <- match.arg(type, type)
  tstat <- testStatUM(y, type=type)
  n <- length(y)
  if (type == "p") {
    Q01 <- -3.110 - 4.652/n - 51.466/(n^2)
    Q05<-  -2.531 - 2.062/n - 17.529/(n^2)
    Q10<-  -2.233 - 1.219/n - 8.178/(n^2)
    }
  else {
    Q01<- -19.51 + 70.09/n - 128.74/(n^2)
    Q05<- -13.02 + 25.05/n - 18.38/(n^2)
    Q10<- -10.1974 + 11.8840/n + 0.8754/(n^2)
    }
  ans<-matrix(c(tstat,Q01,Q05,Q10), nrow=4, ncol=1)
  dimnames(ans)[[1]]<-c("test statistic", 
    "1% critical point", 
    "5% critical point",
    "10% critical point")
  dimnames(ans)[[2]]<-" "
  ans
  }
