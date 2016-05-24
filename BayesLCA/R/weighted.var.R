weighted.var <-
function(x,w,...)	sum(w*(x-weighted.mean(x,w,...))^2)/(sum(w)-1)
