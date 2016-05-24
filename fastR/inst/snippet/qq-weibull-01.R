life01 <- c(16, 34, 53, 75, 93, 120, 150, 191, 240, 339)
qweib <- function(x) { qweibull(x,1.2,146) } 
myplot <- xqqmath(~life01,distribution=qweib,idline=TRUE,qqmathline=FALSE)
