`getJointInf` <-
function(S.matrix, S.observed)
{
    pvalue = 0;
    pvalue.se = 0;
    
    ## GET DISTRIBUTION TABLE
        
    distTable = getDistTable(S.matrix);
    colnames(distTable) = c(names(S.observed),"freq");     
    distribution = distTable;      
    mc.size = nrow(S.matrix);
    
    ## GET P-VALUE ESTIMATE
            
    index = row.matches(S.observed,distTable[,1:length(S.observed)]);
    
    if(length(index) == 0)
    {   
        pvalue = 0;
        pvalue.se = 0;
        
        warning("observed value of the sufficient statistics for 'joint' was not sampled", call. = FALSE);
    }
    else
    {
        observed.freq = distTable[index,ncol(distTable)];
        pvalue = sum(distTable[distTable[,ncol(distTable)] <= observed.freq,ncol(distTable)]);
        
        joint.unique = distTable[,1:length(S.observed)];
        joint.unique.freq = distTable[,ncol(distTable)];
        
        values = rep(0,nrow(S.matrix));
                
        for(j in 1:nrow(joint.unique))
        {
            temp = row.matches(as.vector(joint.unique[j,],mode='numeric'),S.matrix);
            
            values[temp] = joint.unique.freq[j];
        }
        
        values=as.vector(values);
        
        if(nrow(joint.unique) > 1)
        {
            batch.size=batchSize(values);
        
            if(batch.size > length(values)/10)
            {
                batch.size=sqrt(length(values));
            }
            
            pvalue.se = sqrt(bm(vals=values,bs=batch.size,g=function(x) return(as.vector(x <= observed.freq,mode='numeric')))$var);
        }
        else
        {   
            warning("conditional distribution of the joint sufficient statistics was found to be degenerate", call. = FALSE);
                    
            pvalue.se = NA;
        }
    }
    
    joint = list(pvalue,pvalue.se,nrow(S.matrix),distTable);
    
    names(joint) = c("pvalue","pvalue.se","mc.size","distribution");
    
    return(joint);
}

