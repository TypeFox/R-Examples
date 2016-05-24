#get list of intruders (tips) only displayed out of the outputlist of AssessMonophyly
# written by Orlando Schwery 2015

GetIntruderTips <-
function(solution, taxa=NULL, taxlevels='ALL') {
    alltips <- list()  # create empty list to be filled
    if (taxlevels == 'ALL') {  # if all taxlevels are looked for
        for (i in 1:length(solution)) {  # loop through all taxlevels
            nametip <- names(solution)[i]  # create namelabel for current taxlevel
            if (length(taxa) == 0) {  # pick all if no taxon specified
                tmp <- solution[[i]]$IntruderTips  # extract sub-list of intruder tips from solution
                alltips[[nametip]] <- tmp  # add extracted list as sub-object to output list and label it with the appropriate taxlevel nr.
            } else {  # display specific taxon if requested
                alltips2 <- list()  # create empty list to be filled
                for (itx in 1:length(taxa)) {  # loop through vector of taxon names
                    nametip2 <- taxa[itx]  # create label with name of invaded taxon
                    tmp <- solution[[i]]$IntruderTips[[taxa[itx]]]  # extract invaders for this taxon from respective sub-object of solution
                    alltips2[[nametip2]] <- tmp  # add extracted invaders to output list, named after invaded taxon
                }
                alltips[[nametip]] <- alltips2  # add compiled invaders of this taxlevel as sub-list to outputlist, named after current taxlevel
            }
        }
    } # if a specific taxlevel should be focused on
    if (taxlevels != 'ALL' & class(taxlevels) != 'numeric') {  # if named taxlevel requested
      nametip <- paste(taxlevels)  # create namelabel for current taxlevel
      if (length(taxa) == 0){  # pick all if no taxa specified
          tmp <- solution[[taxlevels]]$IntruderTips   # extract sub-list of intruder tips from solution
          alltips[[nametip]] <- tmp  # add extracted list as sub-object to output list and label it with the appropriate taxlevel nr.
      } else {  # display specific taxa if specified
          alltips2 <- list()  # create empty list to be filled
          for (itx in 1:length(taxa)) {  # loop through vector of taxon names
              nametip2 <- taxa[itx]  # create label with name of invaded taxon first
              tmp <- solution[[taxlevels]]$IntruderTips[[taxa[itx]]]  # extract invaders for this taxon from respective sub-object of solution
              alltips2[[nametip2]] <- tmp  # add extracted invaders to output list, named after invaded taxon
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
            tmp <- solution[[taxlevels]]$IntruderTips   # extract sub-list of intruder tips from solution
            alltips[[nametip]] <- tmp  # add extracted list as sub-object to output list and label it with the appropriate taxlevel nr.
        } else {  # display specific taxa if specified
            alltips2 <- list()  # create empty list to be filled
            for (itx in 1:length(taxa)) {  # loop through vector of taxon names
                nametip2 <- taxa[itx]  # create label with name of invaded taxon first
                tmp <- solution[[taxlevels]]$IntruderTips[[taxa[itx]]]  # extract invaders for this taxon from respective sub-object of solution
                alltips2[[nametip2]] <- tmp  # add extracted invaders to output list, named after invaded taxon
            }
            alltips[[nametip]] <- alltips2  # add compiled invaders of this taxlevel as sub-list to outputlist, named after current taxlevel
        }
    }
    if (length(alltips) == 0) {
      stop('No results found. Consider reviewing taxlevels argument.')
    }
    alltips  # export output list
}
