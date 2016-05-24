loadALL <-
function(direc=getwd())
{
if(is.list(direc))
  {
  getitall <- lapply(direc,function(x)
    {
    
    allf <- dir(x) 
    rsourcefiles <- grep("^.*\\.R$",allf,value=TRUE,perl=TRUE)
    rsourcefiles2 <- paste(x,"/",rsourcefiles,sep="")
    }) 
  
  for(i in unlist(getitall))
  {
    source(i)  
  }
  
  } else {
          
          rsourcefiles <- grep("^.*\\.R$",dir(direc),value=TRUE,perl=TRUE)
          rsourcefiles2 <- paste(direc,"/",rsourcefiles,sep="")
          for(i in rsourcefiles2)
            {
            source(i)  
            }
              
          }



}
