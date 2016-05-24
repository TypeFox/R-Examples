#### function for printing the prototest object
#### input:
####    - x = Object of class prototest
print.prototest <- function(x, ...){
  out = data.frame (ts=round(x$ts, 3), p.val=round(x$p.val, 4))
  
  print (out, ...)
}