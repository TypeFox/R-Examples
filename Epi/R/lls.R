lls <-
# A function that expands the functionality of ls()
function( pos = 1, pat = "", all=FALSE, print=TRUE )
{
# First a function that returns length/dim when you ask for it
dimx <- function(dd) if (is.null(dim(dd))) length(dd) else dim(dd)
# A vector of object names
lll <- ls( pos=pos, pattern=pat, all.names=all )
# Are there any objects at all?
if( length(lll) > 0 )
{
obj.mode <-
obj.clas <-
obj.size <- character(0)
# Then find mode, class, name and dimension of them and return it
for(i in 1:length(lll))
{
obj.mode[i] <-        eval( parse(text = paste( "mode(", lll[i], ")")))
obj.clas[i] <- paste( eval( parse(text = paste("class(", lll[i], ")"))), collapse=" " )
obj.size[i] <- paste( eval( parse(text = paste( "dimx(", lll[i], ")"))), collapse=" " )
}
dfr <-
data.frame( name=lll,
            mode=obj.mode,
           class=obj.clas,
            size=obj.size,
stringsAsFactors=FALSE )
print( invisible( dfr ), right=FALSE )
}
}

