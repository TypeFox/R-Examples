 
euclidean.distance.f <- function(fossite, modsite){ 
 
#fossite is a two element vector containing the x and y coordinate of the fossil site
#modsite is a two element vector containing the x and y coordinate of the modern site
# Steven Mosher: added return()
       return( sqrt(sum((fossite - modsite)*(fossite - modsite))))
  }



