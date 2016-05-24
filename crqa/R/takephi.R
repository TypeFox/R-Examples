## extract the phi coefficient from two *already* lagged series
## Arguments: x,y (the two lagged series)
##            k = object of interest

.packageName <- 'crqa'

takephi <- function(x, y, k){
    
    sequence = seq(1,length(x),1)
    CumCT = matrix(0, ncol=2, nrow=2)
    
    ind = grep(k, x); ind.not = setdiff(sequence,ind)
    
    ind1  = ind2 = 0;
    
    for (t in 1:length(ind)){
        
        id = length(grep(paste(k, "$", sep=""), y[ind[t]]))
        
        if (id >= 1){    
            ind1 = ind1 + 1
        } else {
            ind2 = ind2 + 1}
      
    }
    
    CumCT[1,1] = ind1; CumCT[1,2] = ind2;
    
    ind1  = ind2 = 0;
    for (t in 1:length(ind.not)){
        
        id = length(grep(paste(k, "$", sep=""), y[ind.not[t]]))
        
        if (id >= 1){    
            ind1 = ind1 + 1
        } else {
            ind2 = ind2 + 1}
            
    }
    
    CumCT[2,1]=ind1; CumCT[2,2]=ind2;
    
    marg = c(sum(CumCT[1,]),sum(CumCT[2,]),sum(CumCT[,1]),sum(CumCT[,2]));
    
    phi = (CumCT[1,1]*CumCT[2,2]- CumCT[1,2]*CumCT[2,1])/sqrt(prod(marg))
    
    return(phi)
}
