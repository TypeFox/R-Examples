coking <- data.frame(
  width = gl(3,6,18,labels=c(4,8,12)),
  temp  = gl(2,3,18,labels=c(1600,1900)),
  time  = scan("coking.txt",quiet=TRUE)
)