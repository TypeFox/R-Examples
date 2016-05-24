NArray <-
function( x, cells=NA )
{
if( !is.list(x) ) stop("Argument must be a (named) list." )
array( cells, dimnames=x, dim=sapply( x, length ) )
}

ZArray <- function( x, cells=0 ) NArray( x, cells=cells )

larray <-
function( data=NA, dim, dimnames )
{
if( is.list(data) )
  return( array(data=NA,dim=sapply(data,length),
                        dimnames=data) )
else
if( is.list(dim) )
  return( array(data=data,dim=sapply(dim,length),
                          dimnames=dim) )
else
if( is.list(dimnames) )
  return( array(data=data,dim=sapply(dimnames,length),
                          dimnames=dimnames) )
else
if( !missing(dimnames) )
  return( array(data=data,dim=length(data),
                          dimnames=dimnames) )
else
  return( array(data=data,dim=length(data)) )
}
