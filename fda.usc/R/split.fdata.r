#################################################
c.ldata<-function (x, f, drop = FALSE, ...) 
{
    if (!is.list(x)) stop("x is not a list") 
    lenl<-length(x)
    out<-x
     lenf<-length(f)
    for (i in 1:lenl) {
        if (lenf>nrow(x[[i]])) stop("Incorrect length of f")
       out[[i]] <- x[[i]][f,]
    }
    out
}
#################################################
split.fdata<-function(x,f,drop=FALSE,...){
if (!is.factor(f)) f<-factor(f)
nlev<-nlevels(f)
lev<-levels(f)
if (is.matrix(x$data)) x$data<-data.frame(x$data)
if (is.fdata(x)) {
    out<-split(x$data,f,drop=drop,...)
    }
for (i in 1:nlev) out[[lev[i]]]<-fdata(out[[lev[i]]],x$argvals,x$rangeval,x$names)
out
}

#################################################
unlist.fdata<-function(x, recursive = TRUE, use.names = TRUE){
nlev<-length(x)
lev<-names(x)
dat<-x[[1]]$data
arg<-x[[1]]$argvals
rng<-x[[1]]$rangeval     
for (i in 2:nlev){
   dat<-rbind(dat,x[[i]]$data)
   arg<-rbind(arg,x[[i]]$argvals)   
   rng<-rbind(rng,x[[i]]$rangeval)      
}
if (any(diff(arg))>0 ) stop("Error in argvals")
if (any(diff(rng))>0 ) stop("Error in rangeval")
dat<-fdata(dat,arg[1,],rng[1,],x[[1]]$names)
return(dat)
}

