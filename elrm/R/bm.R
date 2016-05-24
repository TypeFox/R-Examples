`bm` <-
function(vals,bs="sqroot",g=id,warn=FALSE)
{
    N <- length(vals);
   
    if(is.null(N))
    {
        if(warn)
        {
            warning("vals must be a vector", call. = TRUE);
        }
    }
    
    if (N<1000)
    {
        if (warn) # if warning
        {
            warning("too few samples (less than 1000)", call. = TRUE);
        }
        if (N<10)
        {
            return(NA);
        }
    }

    if (bs=="sqroot") 
    {
        b <- floor(sqrt(N)); # batch size
        a <- floor(N/b); # number of batches
    }
    else if (bs=="cuberoot") 
    {
        b <- floor(N^(1/3)); # batch size
        a <- floor(N/b); # number of batches
    }
    else # batch size provided
    {
        stopifnot(is.numeric(bs))  
        b <- floor(bs) # batch size
        if (b > 1) # batch size valid
        {
            a <- floor(N/b); # number of batches
        }
        else
        {
            stop("batch size invalid (bs=",bs,")");
        }
    }
    
    if(a < 10)
    {
        stop("too few batches (less than 10)\n");
    }
    
    Ys = sapply(1:a,function(k) {return(mean(g(vals[((k-1)*b+1):(k*b)])))});
    
    muhat <- mean(Ys);
    sigmahat <- var(Ys);
    bmvar <- sigmahat/a;

    return(list(est=muhat,var=bmvar,bs=bs,Ys=Ys));
}

