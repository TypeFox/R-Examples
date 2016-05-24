################################
# Similarity class and methods
################################
# deprecated
## similarity <- function(start1, end1, start2, end2, col=NULL){
##   if (missing(start1) | missing(end1) | missing(start2) | missing(end2))
##     stop("Both starts and ends must be provided")
##   similarity <- list(start1=start1, end1=end1, start2=start2, end2=end2,
##                      col=col)
##   warning("Similarity is deprecated. Use 1-line comparison instead")
##   as.similarity(similarity)
## }
## as.similarity <- function(similarity){
##   # return self is already of the right class
##   if (is.similarity(similarity)) {
##     return(similarity)
##   }
##   # check for presence of arguments
##   if (any(is.null(c(similarity$start1, similarity$end1,
##                     similarity$start2, similarity$end2))))
##     stop("One start or end missing")
##   if (is.null(similarity$col)) similarity$col <- "grey"
##   # check for correct argument types
##   if (any(!is.numeric(c(similarity$start1, similarity$end1,
##                         similarity$start2, similarity$end2))))
##     stop("One start or end missing")
##   # class attribution, overriding precedent classes
##   class(similarity) <- c("similarity", "list")
##   warning("Similarity is deprecated. Use 1-line comparison instead")
##   similarity
## }
## is.similarity <- function(similarity){
##   warning("Similarity is deprecated. Use 1-line comparison instead")
##   inherits(similarity, "similarity")
## }

