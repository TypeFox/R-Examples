concatenate.regions <- function(object){

n.regions <- length(object@region.names)

return(concatenate_to_whole_genome(object,n.regions))

}
