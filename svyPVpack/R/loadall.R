loadall <-
function(patH=NA)
    {
    
    if(is.na(patH))
          {
          patH <- getwd()
          }
    
    nams <- dir(patH)  
      
    for(i in nams)   
          {
          source(paste0(patH,"/",i)) 
          }  
      
      
    }
