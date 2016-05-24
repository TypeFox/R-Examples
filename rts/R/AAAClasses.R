# Author: Babak Naimi, naimi.b@gmail.com
# Date :  July 2012
# Version 1.0
# Licence GPL v3

setOldClass("xts")
setClass("RasterStackTS",
         representation(raster="RasterStack",
                        time="xts"),
         validity=function(object){
           return (nlayers(object@raster) == length(object@time))
         }  
         )


setClass("RasterBrickTS",
         representation(raster="RasterBrick",
                        time="xts"),
         validity=function(object){
           return (nlayers(object@raster) == length(object@time))
         }  
         )


setClassUnion("RasterStackBrickTS", c("RasterStackTS", "RasterBrickTS"))


setClass ('rts', contains = c('xts') )
