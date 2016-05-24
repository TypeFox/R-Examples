# get summary only displayed out of the outputlist of AssessMonophyly
# written by Orlando Schwery 2015

GetSummaryMonophyly <-
function(solution, taxlevels='ALL') {
    Allsummaries <- list()  # create empty list to be filled
    if (taxlevels == 'ALL') {  # if all taxlevels are looked for
        for (i in 1:length(solution)) {  # loop through all taxlevels
            namesum <- names(solution)[i]  # create namelabel for current taxlevel
            tmp <- (solution[[i]]$summary)  # extract summary object from solution
            Allsummaries[[namesum]] <- tmp  # add extracted table as sub-object to output list and label it with the appropriate taxlevel nr.
        }
    }   # if only a specific taxlevel is requested
    if (taxlevels != 'ALL' & class(taxlevels) != 'numeric') {  # if named taxlevel requested
        namesum <- paste(taxlevels)
        tmp <- (solution[[taxlevels]]$summary)
        Allsummaries[[namesum]] <- tmp
    }
    if (class(taxlevels) == 'numeric') {
        if (taxlevels > length(solution)) {  # test whether requested taxlevel is among available ones and display error if not
            stop('Requested taxonomic level not available (less levels specified as analysis input)!')
        }
        namesum <- names(solution)[taxlevels]  # create namelabel for current taxlevel
        tmp <- (solution[[taxlevels]]$summary)  # extract summary object from solution
        Allsummaries[[namesum]] <- tmp  # add extracted list as sub-object to output list and label it with the appropriate taxlevel nr.
    }
    if (length(Allsummaries) == 0) {
      stop('No results found. Consider reviewing taxlevels argument.')
    }
    Allsummaries  # export output list
}
