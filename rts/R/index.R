# Author: Babak Naimi, naimi.b@gmail.com
# Date :  July 2012
# Version 1.0
# Licence GPL v3



`index.RasterStackBrickTS` <- function(x,...){
  index(x@time,...)
}

`index<-.RasterStackBrickTS` <- function(x,value){
  index(x@time) <- value
}


