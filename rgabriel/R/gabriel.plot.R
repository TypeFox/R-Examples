gabriel.plot <- function(x, f, upper, lower=upper, length=0.1,...){
  if (class(f) != "factor") {
    f  <-  factor(f)
  }  
  input  <-  cbind(f[!sapply(is.na(x), all)], x[!sapply(is.na(x), all)])
  f  <-  factor(input[, 1])
  x  <-  input[, 2]
  meang <- tapply(x,f,mean)
  tem <- barplot(t(meang),ylim=c(0,max(meang)+2*max(upper)),...)
  arrows(tem,meang+upper, tem, meang-lower, angle=90, code=3, length=length, ...)
}