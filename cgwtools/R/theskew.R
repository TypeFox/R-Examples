theskew <-
function (x) {
# skew and kurtosis from  the e1071 package.
x<-as.vector(x)
 sum((x-mean(x))^3)/(length(x)*sd(x)^3)
 }
