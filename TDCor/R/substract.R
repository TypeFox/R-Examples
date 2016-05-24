substract <-
function(u,v)

{u[na.omit(match(v,u))]<-NA

 unique(na.omit(u))}
