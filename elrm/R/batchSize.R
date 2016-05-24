`batchSize` <-
function(vec)
{
    N = length(vec);
    
    b <- floor(N^(1/3)); # batch size
    a <- floor(N/b);     # number of batches
    
    func = function(bs)
    {
        batches = bm(vals=vec,bs=round(bs,0),g=id)$Ys;
        ac = acf(x=batches,lag.max=2,plot=F)$acf[2];
        
        return(abs(ac));
    }
    
    if(a > 10)
    {
        lower = b;
        upper = floor(N/10);
        b = optimize(f=func,lower=lower,upper=upper)$minimum;
    }
    
    return(round(b,0));
}

