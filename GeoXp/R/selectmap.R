`selectmap` <-
function(var1,var2,obs,Xpoly,Ypoly,method="")
{

####################################################
# Sélection d'un point
####################################################

if (method == "point")
{
    diff<-abs(var1 - as.numeric(Xpoly)) * (max(var2) - min(var2)) + abs(var2 - as.numeric(Ypoly)) * (max(var1) - min(var1));
    if(min(diff[diff==min(diff)]/(max(var2)-min(var2))/(max(var1) - min(var1)))<0.01)
    {
      if (length(obs) == length(var1))
        {
         obs[diff==min(diff)] <- !obs[diff==min(diff)]  
        }
        else
        {
         obs[diff==min(diff),] <- !obs[diff==min(diff),]   
        }
    }
    return(obs)
  }

####################################################
# Sélection d'un polygone
####################################################

 if (method == "poly") 
 {
  polyg<-cbind(unlist(Xpoly),unlist(Ypoly))
  def <- inout(cbind(var1,var2), polyg, bound=TRUE) 
            
  ifelse(length(obs) == length(var1), obs[def] <- !obs[def], obs[def, ] <- !obs[def, ])
  return(obs)
 }
}

