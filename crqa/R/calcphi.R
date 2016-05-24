# GNU License: written by Moreno I. Coco (moreno.cocoi@gmail.com)

### Calculate the object (recurrence) coefficient (match,mismatch) across lags
##  Arguments: (ts1, ts2) = two categorical time-series
##              k = the object from which the phi is calculated
##              lag = a sequence of lags over which we want to calculate recurrence


.packageName <- 'crqa'

calcphi <- function(t1, t2, ws, k){

    t1 = as.vector(as.matrix(t1)) ## make sure ts are vectors
    t2 = as.vector(as.matrix(t2))
    
    phis = vector();
    
    for (i in (-ws-1):-2){
        
        ix = abs(i);
        y = t2[ix:length(t2)];
        x = t1[1:length(y)];
        
        phis = c(phis, takephi(x, y, k) )
        
    }
    
    phis = c(phis,  takephi(t1, t2, k) ) ## this is phi at lag 0    
    
    
    for (i in 2:(ws+1)){
        
        x = t1[i:length(t1)]
        y = t2[1:length(x)]
    
        phis = c(phis, takephi(x, y, k) )

    }

return( phis )

}


