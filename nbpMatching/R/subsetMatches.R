#'Subset Matches
#'
#'Remove unpaired or unnecessary matches.
#'
#'Given a nonbimatch object, remove elements matched to phantoms, chameleons, or
#'ghosts.  Also remove pairs whose distance is infinite.
#'
#'@aliases subsetMatches subsetMatches,nonbimatch-method
#'@param matches A nonbimatch object.
#'@param phantom A logical value.  Remove elements matched to phantom elements.
#'@param chameleon A logical value.  Remove elements matched to chameleon
#'elements.
#'@param ghost A logical value.  Remove elements matched to ghost elements.
#'@param infinite A logical value.  Remove elements matched at infinite
#'distance. This will include elements forced to match in spite of having an
#'infinite distance set by the prevent option in \code{\link{gendistance}}.
#'@param halvesOnly A logical value.  Use halves element instead of matches.
#'@return a data.frame
#'@exportMethod subsetMatches
#'@author Cole Beck
#'@examples
#'
#'df <- data.frame(id=LETTERS[1:25], val1=rnorm(25), val2=rnorm(25))
#'df.dist <- gendistance(df, idcol=1, ndiscard=4)
#'df.mdm <- distancematrix(df.dist)
#'df.match <- nonbimatch(df.mdm)
#'subsetMatches(df.match)
#'subsetMatches(df.match, halvesOnly=FALSE)
#'subsetMatches(df.match, phantom=FALSE)
#'

setGeneric("subsetMatches", function(matches, phantom=TRUE, chameleon=TRUE, ghost=TRUE, infinite=TRUE, halvesOnly=TRUE) {
    standardGeneric("subsetMatches")
})

setMethod("subsetMatches", signature(matches="nonbimatch"), function(matches, phantom=TRUE, chameleon=TRUE, ghost=TRUE, infinite=TRUE, halvesOnly=TRUE) {
    if(halvesOnly) {
        dat <- matches$halves
    } else {
        dat <- matches$matches
    }
    id1 <- dat[,'Group1.ID']
    id2 <- dat[,'Group2.ID']
    ix <- numeric()
    if(phantom) {
        ix <- union(ix, c(grep('^phantom', id1), grep('^phantom', id2)))
    }
    if(chameleon) {
        ix <- union(ix, c(grep('^chameleon', id1), grep('^chameleon', id2)))
    }
    if(ghost) {
        ix <- union(ix, c(grep('^ghost', id1), grep('^ghost', id2)))
    }
    if(infinite) ix <- union(ix, which(is.infinite(dat[,'Distance'])))
    if(length(ix)) dat <- dat[-ix,]
    dat
})
