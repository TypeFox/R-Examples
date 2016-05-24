thekurt <-
function (x) {
# skew and kurtosis from  the e1071 package.
x<-as.vector(x)
  sum((x-mean(x))^4)/(length(x)*var(x)^2) - 3
  }
