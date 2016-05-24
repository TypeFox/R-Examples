#' Make ztable from object mytable
#'
#'@param x An object of mytable
#'@param digits Numeric vector of length equal to one (in which case it will be
#'       replicated as necessary) or to the number of columns of the resulting table
#' @param ... arguments to be passed to \code{\link{ztable_sub}}
#'@examples
#'require(moonBook)
#'res=mytable(sex~.,data=acs)
#'z=ztable(res)
#'z
ztable.mytable=function(x,digits=NULL,...){

    count=ncol(x$res)
    if(x$show.all==FALSE) count=ncol(x$res)-7
    myalign="ll"
    for(i in 2: count) myalign=paste(myalign,"c",sep="")
    z=ztable(x$res[1:count],align=myalign)
    colnames(z$x)[1]=""
    sub=paste("(N=",x$count,")",sep="")
    sub=c("",sub)
    while(length(sub)<count) sub=c(sub,NA)
    z=addSubColNames(z,sub)
    z$include.rownames=FALSE
    z=vlines(z,type=0)
    z
}
