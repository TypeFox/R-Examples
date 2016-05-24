p2chi = function(p, k, steps = 30, chimax = 1E31){
    
    temp = function(p, k, steps, chimax){
        
        chiseq = seq(0, chimax, len=100)
        pseq = chipval(X=chiseq, k=k)
        chimid = which.min(abs(p-pseq))
        
        for(i in 1:steps){
            
            chimin = max((chimid-1), 1)
            chimax = min((chimid+1), 100)
            chiseq = seq(chiseq[chimin], chiseq[chimax], len=100)
            pseq = chipval(X=chiseq, k=k)
            chimid = which.min(abs(p-pseq))
            
        }
        
        return(chiseq[chimid])
        
    }
    
    return(unlist(lapply(p, temp, k=k, steps=steps, chimax=chimax)))
    
}
