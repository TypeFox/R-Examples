#' Make ztable from object cbind.mytable
#'
#'@param x An object of cbind.mytable
#'@param digits Numeric vector of length equal to one (in which case it will be
#'       replicated as necessary) or to the number of columns of the resulting table
#'@param ... arguments to be passed to \code{\link{ztable_sub}}
#'@examples
#'require(moonBook)
#'res=mytable(sex+DM~.,data=acs)
#'z=ztable(res)
#'z
ztable.cbind.mytable=function(x,digits=NULL,...){
    t=list()
    myalign="ll"
    for(i in 1:length(x)){
        count=ncol(x[[i]]$res)
        if(x[[i]]$show.all==FALSE) count=ncol(x[[i]]$res)-7

        for(j in 2: count) myalign=paste(myalign,"c",sep="")
        if(i==1) t[[i]]=x[[i]]$res[1:count]
        else t[[i]]=x[[i]]$res[2:count]
        sub1=paste("(N=",x[[i]]$count,")",sep="")
        if(i==1) {
            sub1=c(NA,sub1)
            while(length(sub1)<count) sub1=c(sub1,NA)
        } else {
            while(length(sub1)<(count-1)) sub1=c(sub1,NA)
        }

        if(i==1) sub=sub1
        else sub=c(sub,sub1)
    }
    mydf=t[[1]]
    for(i in 2:length(x)) mydf=cbind(mydf,t[[i]])
    caption=paste("Descriptive Statistics Stratified by \'",
                  toupper(attr(x,"group")[1]),"\' and \'",
                  toupper(attr(x,"group")[2]),"\'",sep="")
    z=ztable(mydf,caption=caption,align=myalign)
    z=addSubColNames(z,sub)
    z$include.rownames=FALSE
    #colnames(z$x)[1]=""
    #cgroup=c(toupper(attr(x,"group")[1]),attr(x,"caption"))
    cgroup=c("",attr(x,"caption"))
    colnames(z$x)[1]=""
    n.cgroup=c(1,rep(count-1,length(x)))
    z=addcgroup(z,cgroup=cgroup,n.cgroup)
    z=vlines(z,type=0)
    z
}


