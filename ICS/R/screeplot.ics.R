### screeplot method for an ics object
### barplot or lineplot are options
###

`screeplot.ics` <-
    function(x,index=NULL,type="barplot",
               main = deparse(substitute(x)),ylab="generalized kurtosis",xlab= "component",...)
    {
    if(class(x)!="ics") stop("'x' must be of class ics")
    type<-match.arg(type,c("barplot","lines"))
    if (is.null(index)) index=1:length(x@gKurt)
    if (type=="barplot")
        {
        barplot(x@gKurt[index],names.arg=index,ylab=ylab,xlab=xlab,main=main,...)
        }
    else
        {
        plot(index,x@gKurt[index],type="b",ylab=ylab,xlab=xlab,axes=F,main=main,...)
        axis(2)
        axis(1, at = index, labels = T)
        }
    invisible()
    }
