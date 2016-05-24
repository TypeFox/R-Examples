`allignment` <-
function(response, block, type=c("mean","median"),...)
{

args<-list(...)

if(length(response)!=length(block))
 {stop("'response' and 'block' must be vectors of the same length")}

if(!is.numeric(response))
 {stop("'response' must be a numeric vector")}

type <- match.arg(type)

f <- as.factor(block)

splitdat <- split(response, f=f)

out1<-lapply(X=splitdat,
 FUN=function(x,...)
{
args$x<-x
m<-do.call(type, args)
return(x-m)
}
)

out<-unsplit(value=out1, f=f)
return(out)

}

