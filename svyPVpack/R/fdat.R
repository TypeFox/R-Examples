
fALL <- function(mybys,pm,svydat)
{
  
  

facordall <- mapply(function(x,number){
  
  if(is.factor(pm[,number]))
  {
    #factor(pm[,number], levels=levels(svydat$variables[[x]]))
    factor(pm[,number], levels=levels(pm[,number]))
  } else 
  {
    #factor(pm[,number], levels=levels(as.factor(svydat$variables[[x]])))
    factor(pm[,number], levels=levels(as.factor(pm[,number])))
  }
}, x=mybys, number=1:length(mybys),SIMPLIFY=FALSE)


return(data.frame(facordall))
  
}











