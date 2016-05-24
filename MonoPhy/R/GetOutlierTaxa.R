#get list of outliers (genera) only displayed out of the outputlist of AssessMonophyly
# written by Orlando Schwery 2015

GetOutlierTaxa <-
function(solution, taxlevels='ALL') {
    alltaxa <- list()  # create empty list to be filled
    if (taxlevels == 'ALL') {  # if all taxlevels are looked for
        for (i in 1:length(solution)) {  # loop through all taxlevels
            nametax <- names(solution)[i]  # create namelabel for current taxlevel
            tmp <- solution[[i]]$OutlierTaxa  # extract sub-list of outlier taxa from solution
            alltaxa[[nametax]] <- tmp # add extracted list as sub-object to output list and label it with the appropriate taxlevel nr.
        }
    }  # if a specific taxlevel should be focused on
    if (taxlevels != 'ALL' & class(taxlevels) != 'numeric') {  # if named taxlevel requested
      nametax <- paste(taxlevels)  # create namelabel for current taxlevel
      tmp <- solution[[taxlevels]]$OutlierTaxa  # extract sub-list of outlier taxa from solution
      alltaxa[[nametax]] <- tmp  # add extracted list as sub-object to output list and label it with the appropriate taxlevel nr.
    }
      if (class(taxlevels) == 'numeric') {
        if (taxlevels > length(solution)) {  # test whether requested taxlevel is among available ones and display error if not
            stop('Requested taxonomic level not available (less levels specified as analysis input)!')
        }
		    nametax <- names(solution)[taxlevels]  # create namelabel for current taxlevel
        tmp <- solution[[taxlevels]]$OutlierTaxa  # extract sub-list of outlier taxa from solution
        alltaxa[[nametax]] <- tmp  # add extracted list as sub-object to output list and label it with the appropriate taxlevel nr.
    }
    if (length(alltaxa) == 0) {
      stop('No results found. Consider reviewing taxlevels argument.')
    }
    alltaxa  # export output list
}
