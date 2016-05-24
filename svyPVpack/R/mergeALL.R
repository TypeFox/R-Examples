mergeALL <-
function(INP)
{

rootf <- INP[[1]]   
  
for(i in 2:length(INP))  
    {
      
    rootf <- merge(rootf,INP[[i]],all=TRUE, sort=FALSE) 

    }

rootf  
}
