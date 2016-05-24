MS_getStats <- function(object,locus=1,population=1){
# object is an object of class cs.stats
return(object@locus[[locus]]@stats[[population]])
}
