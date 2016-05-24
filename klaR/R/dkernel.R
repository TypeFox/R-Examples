dkernel<-function(x, kernel=density(x), interpolate=FALSE, ...)
{
foo<-function(x,kernel,n)
{
if (x <= kernel$x[1]) return(kernel$y[1])
else if (x >= kernel$x[n]) return(kernel$y[n])
else 
    {
        pos<-which.min((kernel$x-x)^2)
        if (kernel$x[pos]>x) return(mean(c(kernel$y[pos],kernel$y[(pos-1)])))
        else return(mean(c(kernel$y[pos],kernel$y[(pos+1)])))
    }
}
n<-length(kernel$x)
if (interpolate) y<-sapply(x,foo,kernel=kernel,n=n)
else y<-sapply(x,FUN=function(y){kernel$y[(which.min((kernel$x-y)^2))]})
return(y)
}
