`plot.elrm` <-
function(x,p=1.0,breaks="Sturges", ask=FALSE, ...)
{   
    if(p <= 0 || p > 1)
    {
        stop("'","p","'"," must be a number greater than 0 and smaller than or equal to 1. \n");
    }
    
    tmp.ask = options("device.ask.default")[[1]];
    options("device.ask.default"=ask);
    
    S = as.data.frame(x$mc);
    
    n = p*nrow(S);
    
    rsample = sort(sample(nrow(S),n));

    if(!is.character(breaks))
    {
        if(length(breaks) != length(names(S)))
        {
            breaks = rep(breaks[1],length(names(S)));
        }
    } 
    
    for(i in 1:length(names(S)))
    {
        dev.new();
		
        par(mfrow=c(2,1));
        
        plot(y=S[rsample,i],x=rsample,col=(i%%2)+3,pch=19,ylab=as.list(names(S))[i],xlab="iterations",main=paste("Trace for ", names(S)[i]));
        
        if(is.character(breaks))
        {
            hist(S[rsample,i],xlab=names(S[i]),ylab="counts",col=(i%%2)+3,main=paste("Histogram for ",names(S)[i]),breaks=breaks);
        }
        else
        {   
            hist(S[rsample,i],xlab=names(S[i]),ylab="counts",col=(i%%2)+3,main=paste("Histogram for ",names(S)[i]),breaks=breaks[i]);
        }
    }
    
    options("device.ask.default"=tmp.ask);
}

