## Gap Statistics
CalculateGapK <- function(W, WStar, kG, n.monterCarlo)
{
    WLog <- log(W)
    WStarLog <- log(WStar)
    Gap <- apply(WStarLog, 2, mean) - WLog
    s <- apply( WStarLog,2,sd) * sqrt( (n.monterCarlo-1)/n.monterCarlo ) * sqrt(1+1/n.monterCarlo)

    I <- 2;
    gap_kHat <- 1
    for (k in I:(kG-1))
    {
        if ( Gap[k] >= Gap[k+1]-s[k+1] ) 
    	{
    	    gap_kHat <- k
    	    break
    	}
    }

    gap_kHat
}
