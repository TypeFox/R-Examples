AICpls <- function(ncomp,residpls,weights=rep.int(1,length(residpls))) {
if(is.null(weights)){rep.int(1,length(residpls))}
return(-2*loglikpls(residpls,weights)+2*(ncomp+1+1))}
