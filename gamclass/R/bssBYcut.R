bssBYcut <-
function(x, y, data){
xnam <- deparse(substitute(x))
ynam <- deparse(substitute(y))
xv <- data[,xnam]
yv <- data[,ynam]
sumss <- function(x, y, cut){
  av <- mean(y)
  left <- x<cut
  sum(left)*(mean(y[left])-av)^2+sum(!left)*(mean(y[!left])-av)^2
}
xOrd <- unique(sort(xv))[-1]
bss <- numeric(length(xOrd))
for(i in 1:length(xOrd)){
  bss[i] <- sumss(xv, yv, xOrd[i])
  }
data.frame(xOrd=xOrd, bss=bss)
 }
