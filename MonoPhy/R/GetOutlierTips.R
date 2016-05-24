#get list of outliers (tips) only displayed out of the outputlist of AssessMonophyly
# written by Orlando Schwery 2015

GetOutlierTips <-
function(solution, taxa=NULL, taxlevels='ALL') {
    alltips <- list()  # create empty list to be filled
    if (taxlevels == 'ALL') {  # if all taxlevels are looked for
        for (i in 1:length(solution)){  # loop through all taxlevels
            nametip <- names(solution)[i]  # create namelabel for current taxlevel
            if (length(taxa) == 0) {  # pick all if no taxon specified
                tmp <- solution[[i]]$OutlierTips  # extract sub-list of outlier tips from solution
                alltips[[nametip]] <- tmp  # add extracted list as sub-object to output list and label it with the appropriate taxlevel nr.
            } else {  # display specific taxon if specified
                alltips2 <- list()  # create empty list to be filled
                for (i in 1:length(taxa)) {  # loop to go through vector of taxon names
                    nametip2 <- taxa[i]  # create label with name of taxon with outliers first
                    tmp <- solution[[i]]$OutlierTips[[taxa[i]]]  # extract outliers for this taxon from respective sub-object of solution
                    alltips2[[nametip2]] <- tmp  # add extracted outliers to output list, named after the that taxon
                }
                alltips[[nametip]] <- alltips2  # add compiled outliers of this taxlevel as sub-list to outputlist, named after current taxlevel
            }
        }
    }  # if a specific taxlevel should be focused on
    if (taxlevels != 'ALL' & class(taxlevels) != 'numeric') {  # if named taxlevel requested
        nametip <- paste(taxlevels)  # create namelabel for current taxlevel
        if (length(taxa) == 0){  # pick all if no taxa specified
            tmp <- solution[[taxlevels]]$OutlierTips   # extract sub-list of outlier tips from solution
            alltips[[nametip]] <- tmp  # add extracted list as sub-object to output list and label it with the appropriate taxlevel nr.
        } else {  # display specific taxa if requested
          alltips2 <- list()  # create empty list to be filled
          for (i in 1:length(taxa)) {  # loop through vector of taxon names
              nametip2 <- taxa[i]  # create label with name of taxon with outliers first
              tmp <- solution[[taxlevels]]$OutlierTips[[taxa[i]]]  # extract outliers for this taxon from respective sub-object of solution
              alltips2[[nametip2]] <- tmp  # add extracted outliers to output list, named after that taxon
            }
            alltips[[nametip]] <- alltips2  # add compiled invaders of this taxlevel as sub-list to outputlist, named after current taxlevel
        }
      }
      if (class(taxlevels) == 'numeric') {
        if (taxlevels > length(solution)) {  # test whether requested taxlevel is among available ones and display error if not
            stop('Requested taxonomic level not available (less levels specified as analysis input)!')
        }
        nametip <- names(solution)[taxlevels]  # create namelabel for current taxlevel
        if (length(taxa) == 0){  # pick all if no taxa specified
            tmp <- solution[[taxlevels]]$OutlierTips   # extract sub-list of outlier tips from solution
            alltips[[nametip]] <- tmp  # add extracted list as sub-object to output list and label it with the appropriate taxlevel nr.
        } else {  # display specific taxa if requested
            alltips2 <- list()  # create empty list to be filled
            for (i in 1:length(taxa)) {  # loop through vector of taxon names
				        nametip2 <- taxa[i]  # create label with name of taxon with outliers first
                tmp <- solution[[taxlevels]]$OutlierTips[[taxa[i]]]  # extract outliers for this taxon from respective sub-object of solution
                alltips2[[nametip2]] <- tmp  # add extracted outliers to output list, named after that taxon
              }
              alltips[[nametip]] <- alltips2  # add compiled invaders of this taxlevel as sub-list to outputlist, named after current taxlevel
		      }
        }
    if (length(alltips) == 0) {
      stop('No results found. Consider reviewing taxlevels argument.')
    }
    alltips
}
