`m4plSummary`<-
function(X, ...) {
 if (is.data.frame(X))                res <- m4plNoMoreSummary(X)
 if (!is.data.frame(X) && is.list(X)) res <- m4plMoreSummary(  X, ...)
 return(res)
 }

