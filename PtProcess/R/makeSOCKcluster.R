#   temporarily included into PtProcess
#   code becomes broken if one simply changes from snow to parallel
makeSOCKcluster <- function(names, ...)
    parallel::makePSOCKcluster(names, ...)


