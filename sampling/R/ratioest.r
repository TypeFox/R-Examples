ratioest<-function(y,x,Tx,pik)
{if (any(is.na(pik))) 
        stop("there are missing values in pik")
    if (any(is.na(y))) 
        stop("there are missing values in y")
    if (any(is.na(x))) 
        stop("there are missing values in x")
    if (length(y) != length(pik) | length(x)!=length(pik) | length(x)!=length(y)) 
        stop("y, x and pik have different lengths")
    sum(y/pik)*Tx/sum(x/pik)
}


