addSightings <- function (capthist, unmarked = NULL, nonID = NULL, verify = TRUE, ...) {
    if (ms(capthist)) {
        if (is.character(unmarked)) {
            ##pre-read and split text file as list
            unmarked <- read.table(unmarked, ...)
            session <- unmarked[,1]
            unmarked <- split(unmarked[,-1], session)
        }
        if (is.character(nonID)) {
            ##pre-read and split text file as list
            nonID <- read.table(nonID, ...)
            session <- nonID[,1]
            nonID <- split(nonID[,-1], session)
        }
        if (!is.null(unmarked)) {
            if (!inherits(unmarked, 'list'))
                stop ("unmarked should be a list for multi-session capthist")
            if (length(unmarked) != length(capthist))
                stop("number of sessions differs between unmarked and capthist")
        }
        if (!is.null(nonID)){
            if (!inherits(unmarked, 'list'))
                stop ("nonID should be a list for multi-session capthist")
            if (length(nonID) != length(capthist))
                stop("number of sessions differs between nonID and capthist")
        }
        for (i in 1:length(capthist)) {
#             s <- ncol(capthist[[i]])
#             capthist[[i]] <- addSightings(capthist[[i]], unmarked[[i]][,1:s], nonID[[i]][,1:s], 
#                                           verify = FALSE)
            capthist[[i]] <- addSightings(capthist[[i]], unmarked[[i]], nonID[[i]], 
                                          verify = FALSE)
        }
        if (verify) verify(capthist)
        capthist
    }
    else {
        if (is.null(markocc(traps(capthist))))
            stop("traps(capthist) lacks markocc attribute; provide this first")
        S <- ncol(capthist)
        K <- ndetector(traps(capthist))
      
        if (is.character(unmarked))
            ## discard session column as unused
            unmarked <- read.table(unmarked, ...)[1:K,2:(S+1)]
        if (is.character(nonID))
            ## discard session column as unused
            nonID <- read.table(nonID, ...)[1:K,2:(S+1)]
        
        Tu(capthist) <- as.matrix(unmarked)
        Tm(capthist) <- as.matrix(nonID)
        if (verify) verify(capthist)
        capthist
    }
}