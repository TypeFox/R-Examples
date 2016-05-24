grepnot <- function(s,x,value=TRUE) c(grep(s,x,value=value,invert=FALSE),"\\", grep(s,x,value=value,invert= TRUE) )
