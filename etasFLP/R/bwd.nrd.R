bwd.nrd <-
function (x,w=replicate(length(x),1)
,d=2) 
{
    if (length(x) < 2L) 
        stop("need at least 2 data points")
     m<-           weighted.mean(x,w)
     return(sqrt(weighted.mean((x-m)^2,w)) * (length(x)*(d+2)/4)^(-1/(d+4)))
}

