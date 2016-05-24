loglikpls <- function(residpls,weights=rep.int(1,length(residpls))) {
if(is.null(weights)){rep.int(1,length(residpls))}
nn=length(residpls)
return(1/2*(sum(log(weights))-nn*(log(2*pi) + 1 - log(nn) + log(sum(weights*residpls^2)))))}
