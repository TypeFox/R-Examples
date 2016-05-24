rlogit <-
function(x) return(ifelse(x>100,1,exp(x)/(1+exp(x))))
