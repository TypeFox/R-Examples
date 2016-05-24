 
euclidean.compy.f <- function(modsite, fossite){ 
 
#modsite is a two element vector containing the x and y coordinate for the modern pollen site
#fossite is a two element vector containing the x and y coordinate for the fossil pollen site
# Steve Mosher : error in orginal code, returned the x variable difference
        #fossite[1] - modsite[1] 
     return(fossite[2]- modsite[2])
        #get the x-component of the vector (length in x in xy units)
  }



