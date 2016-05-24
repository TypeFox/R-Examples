# get result table only displayed out of the outputlist of AssessMonophyly
# written by Orlando Schwery 2015

GetResultMonophyly <-
function(solution, taxlevels='ALL') {
    Allresults <- list()  # create empty list to be filled
    if (taxlevels == 'ALL') {  # if all taxlevels are looked for
        for (i in 1:length(solution)) {  # loop through all taxlevels
            nameres <- names(solution)[i]  # create namelabel for current taxlevel
            tmp <- (solution[[i]]$result)  # extract result object from solution
            Allresults[[nameres]] <- tmp  # add extracted table as sub-object to output list and label it with the appropriate taxlevel nr.
        }
    }  # if only a specific taxlevel is requested
    if (taxlevels != 'ALL' & class(taxlevels) != 'numeric') {  # if named taxlevel requested
        nameres <- paste(taxlevels)
        tmp <- (solution[[taxlevels]]$result)
        Allresults[[nameres]] <- tmp
    }
    if (class(taxlevels) == 'numeric') {
        if (taxlevels > length(solution)) {  # test whether requested taxlevel is among available ones and display error if not
            stop('Requested taxonomic level not available (less levels specified as analysis input)!')
        }
        nameres <- names(solution)[taxlevels]  # create namelabel for current taxlevel
        tmp <- (solution[[taxlevels]]$result)  # extract result object from solution
        Allresults[[nameres]] <- tmp  # add extracted list as sub-object to output list and label it with the appropriate taxlevel nr.
    }
    if (length(Allresults) == 0) {
      stop('No results found. Consider reviewing taxlevels argument.')
    }
    Allresults  # export output list
}
