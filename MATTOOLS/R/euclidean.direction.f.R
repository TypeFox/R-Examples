 
euclidean.direction.f <- function(modsite, fossite){ 
 
#modsite is a two element vector containing the x and y coordinate for the modern pollen site
#fossite is a two element vector containing the x and y coordinate for the fossil pollen site
# clean up adding return() Steven Mosher change to atan2
        x <- fossite[1] - modsite[1]    # get the x-component of the vector
        y <- fossite[2] - modsite[2]    # get the y-component of the vector
        return(((atan2(y, x) * 180)/pi) - 90)   
        # get the quadrat correct arctangent of the vector and return this
  }



