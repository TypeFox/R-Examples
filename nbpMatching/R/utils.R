# extra utility functions

# list of names for created/fake elements found in matched dataset
created.names <- c("phantom", "ghost", "chameleon")

#'Get named sets of matches
#'
#'Create a factor variable using the names from a matched data set.
#'
#'Calculate a name for each pair by using the ID columns from the matched data
#'set.  Return a factor of these named pairs.
#'
#'@aliases get.sets get.sets,data.frame-method get.sets,nonbimatch-method
#'@param matches A data.frame or nonbimatch object.  Contains information on
#'how to match the covariate data set.
#'@param remove.unpaired A boolean value.  The default is to remove elements
#'matched to phantom elements.
#'@param \dots Additional arguments, not used at this time.
#'@return a factor vector
#'@exportMethod get.sets
#'@author Jake Bowers, \url{http://jakebowers.org/}, Cole Beck
#'@examples
#'
#'df <- data.frame(id=LETTERS[1:25], val1=rnorm(25), val2=rnorm(25))
#'df.dist <- gendistance(df, idcol=1)
#'df.mdm <- distancematrix(df.dist)
#'df.match <- nonbimatch(df.mdm)
#'get.sets(df.match)
#'get.sets(df.match$matches)
#'# include the phantom match
#'get.sets(df.match$matches, FALSE)
#'
setGeneric("get.sets", function(matches, remove.unpaired=TRUE, ...) standardGeneric("get.sets"))
setMethod("get.sets", "data.frame", function(matches, remove.unpaired=TRUE, ...) {
    # thanks to Jake Bowers for providing this function
    sets <- matches[,grep("ID", names(matches))]
    f.sets <- apply(sets, MARGIN=1, FUN=function(x) paste(sort(x), collapse='-'))
    names(f.sets) <- sets[,1]
    if(remove.unpaired) f.sets <- f.sets[grep(paste(created.names, collapse="|"), f.sets, invert=TRUE)]
    factor(f.sets)
})

setMethod("get.sets", "nonbimatch", function(matches, remove.unpaired=TRUE, ...) {
    get.sets(matches$matches, remove.unpaired, ...)
})

#'Calculate scalar distance
#'
#'Calculate the scalar distance between elements of a matrix.
#'
#'Take the absolute difference between all elements in a vector, and return a
#'matrix of the distances.
#'
#'@aliases scalar.dist scalar.dist,vector-method
#'@param x A vector of numeric values.
#'@param \dots Additional arguments, not used at this time.
#'@return a matrix object
#'@exportMethod scalar.dist
#'@author Jake Bowers, \url{http://jakebowers.org/}, Cole Beck
#'@examples
#'
#'scalar.dist(1:10)
#'
setGeneric("scalar.dist", function(x, ...) standardGeneric("scalar.dist"))
setMethod("scalar.dist", "vector", function(x, ...) {
    # thanks to Jake Bowers for providing this function
    if(!is.numeric(x)) stop("x should be numeric")
    outer(x, x, FUN=function(i,j) abs(i-j))
})
