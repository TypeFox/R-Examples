# Little functions to make things work

setDefault<-function(x,val)
  ifelse(is.null(x),val,x)