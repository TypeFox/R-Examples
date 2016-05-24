# This is a frequently used function anyway
se<- function(x,na.rm=TRUE) sqrt(var(x,na.rm=na.rm) /
                            length(x[complete.cases(x)]))
