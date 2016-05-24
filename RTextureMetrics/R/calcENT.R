calcENT <-
function(rawmat)
{
    lnrawmat<-log(rawmat)
    
    size<-dim(lnrawmat)[1]
    
    for(i in 1:size)   #by row
    {
        for(a in 1:size)
        {   
            if(lnrawmat[a,i]=="-Inf")
            {
               lnrawmat[a,i]<-0  
            }                    
        } 
    } 
 
    return(sum(rawmat*lnrawmat)) 
}
