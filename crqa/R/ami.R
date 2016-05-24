#' Average Mutual Information of two time-seroes
#' 
#' @param ts1 first time series
#' @param ts2 second time series.
#' @param lag a vector of the time lags
#' @return a vector of average mutual information. 

#' @examples
#'
#' data(crqa)
#' lag = seq(1,100,1)
#' ts1 = leftmov; ts2 = rightmov
#' ami(ts1, ts2, lag)

#' A peak in v for lag > 0 means ts2 is leading ts1.) 
#' v returns how many bits ts1 and ts2 have in common relative to  
#' the number of bits needed to obtain a binned representation of ts1 or ts2.
#' This is done to make the result close to independent bin size.
#'
#'  Optimal binning is obtained by transforming ts1 and ts2 into percentiles 
#'
#' (Based on ami.m MATLAB function written by Aslak Grinsted)

.packageName <- 'crqa'

ami <- function(ts1, ts2, lag){

    ## a default for lag
    ## lag = 0:min(n/2-1,20);
    ## make sure that lags are integers
    lag = round(lag);

    ts1 = as.vector(ts1)
    ts2 = as.vector(ts2)

    n = length(ts1)

    if (n != length(ts2)){
        stop('ts1 and ts2 should have the same length');
    }

    ts1 = ts1-min(ts1)   
    ts1 = ts1*(1-eps(1))/max(ts1)

    ts2 = ts2-min(ts2)
    ts2 = ts2*(1-eps(1))/max(ts2)

    v = vector(mode = "numeric", length = length(lag))
    
    lastbins = 0

    ii = 2

    for (ii in 1:length(lag)){

        abslag = abs(lag[ii]);
    
        ## Define the number of bins
        bins = floor(1+log2(n-abslag)+0.5)

        if (bins != lastbins){
            bints1 = floor(ts1*bins)+1
            bints2 = floor(ts2*bins)+1
        }
        
        lastbins = bins

        PtS = matrix(0, nrow = bins, ncol = bins)
    
	for (jj in 1:(n-abslag)){
            kk = jj + abslag

            if (lag[ii] < 0) { 
                temp = jj; jj = kk; kk = temp; ## swap 
            }
            
            PtS[bints1[kk], bints2[jj]] = PtS[bints1[kk],bints2[jj]] + 1;
        }
        
        PtS = PtS/(n-abslag)
        PtS = PtS + eps(1) ## avoid division and log of zero

        Pts1 = rowSums(PtS)
        Pts2 = colSums(PtS)
 
        q = PtS/ outer(Pts1,Pts2);    
        q = PtS * log2(q);
    
        v[ii] = sum(q)/log2(bins);

    }

    return(v)
        
}
