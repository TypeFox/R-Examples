`qm4pl` <-
function(p=0.05, S=0, C=0, D=0, s=1/1.702, b=0, c=0, d=1, lower.tail=TRUE, log.p=FALSE) {
 if (log.p      == TRUE)  p <- exp(p)
 if (lower.tail == FALSE) p <-    1-p
 result <- log(((d-D)-(C+c))/(p-(C+c)) - 1)/(-1/sqrt(s^2 + S^2)) + b
 return(result)
 }

