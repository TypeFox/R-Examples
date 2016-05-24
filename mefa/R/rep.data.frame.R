`rep.data.frame` <-
function(x, ...)
as.data.frame(lapply(x, rep, ...))
