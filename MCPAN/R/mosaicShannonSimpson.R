#

mosaicdiv<-function(x, decreasing=NULL, ...)
{
if(!class(x)%in%c("Shannonci","Simpsonci")){warning("This is a convenience function for objects of class 'Shannonci' or 'Simpsonci'!")}

tab<-x$sample.estimate$table

aargs<-list(...)

if(is.null(aargs$las))
 {aargs$las<-1}

if(is.null(aargs$main))
 {aargs$main<-""}


if(is.null(decreasing))
{
aargs$x<-tab
do.call("mosaicplot", aargs)
}
else{

if(!is.logical(decreasing)) 
 {stop("decreasing should be either TRUE, FALSE or NULL")}

if(decreasing)
{
cs<-apply(tab, 2, sum)
CO<-order(cs, decreasing=TRUE)
aargs$x<-tab[,CO]
do.call("mosaicplot", aargs)
}


if(!decreasing)
{
cs<-apply(tab, 2, sum)
CO<-order(cs, decreasing=FALSE)
aargs$x<-tab[,CO]
do.call("mosaicplot", aargs)
}
}

}

