intersect <-
function(u,v)

{unique(u[na.omit(match(v,u))])}
