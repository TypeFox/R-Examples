#get list of MRCA nodes of desired genera only displayed out of the outputlist of AssessMonophyly
# written by Orlando Schwery 2015

GetAncNodes <-
function(solution, taxa=NULL, taxlevels='ALL') {
    allnodes <- list() 	# create empty list to be filled
    if (taxlevels == 'ALL') {  #if no specific taxlevel should be focused on
        for (i in 1:length(solution)){  # for every taxlevel in the present solution
            namenod <- names(solution)[i]  # create namelabel for current taxlevel
            if (length(taxa) == 0){  # display all if no genera specified
                tmp <- solution[[i]]$result[, "MRCA", drop=FALSE]  # extract MRCA column out of result table
                allnodes[[namenod]] <- tmp  # add the extracted table as sub-entry to output list, named according to taxlevel
            } else {  # display specific genera if specified
                tmp <- solution[[i]]$result[taxa, "MRCA", drop=FALSE]  # extract MRCA column out of result table for specified taxa only
                allnodes[[namenod]] <- tmp  # add the extracted table as sub-entry to output list, named according to taxlevel
            }
        }
    }  # if a specific taxlevel should be focused on
    if (taxlevels != 'ALL' & class(taxlevels) != 'numeric') {  # if named taxlevel requested
        namenod <- paste(taxlevels)  # create namelabel for current taxlevel
		    if (length(taxa) == 0){  # display all if no genera specified
          tmp <- solution[[taxlevels]]$result[, "MRCA", drop=FALSE]  # extract MRCA column out of result table
          allnodes[[namenod]] <- tmp  # add the extracted table as sub-entry to output list, named according to taxlevel
        } else {  # display specific genera if specified
          tmp <- solution[[taxlevels]]$result[taxa, "MRCA", drop=FALSE]  # extract MRCA column out of result table for specified taxa only
          allnodes[[namenod]] <- tmp  # add the extracted table as sub-entry to output list, named according to taxlevel
        }
    }
     if (class(taxlevels) == 'numeric') {
        if (taxlevels > length(solution)){  # test if taxlevel number exceeds number of present taxlevels in solution and display error if the case
            stop('Requested taxonomic level not available (less levels specified as analysis input)!')
        }
		    namenod <- names(solution)[taxlevels]  # create namelabel for current taxlevel
		    if (length(taxa) == 0){  # display all if no genera specified
          tmp <- solution[[taxlevels]]$result[, "MRCA", drop=FALSE]  # extract MRCA column out of result table
          allnodes[[namenod]] <- tmp  # add the extracted table as sub-entry to output list, named according to taxlevel
        } else {  # display specific genera if specified
          tmp <- solution[[taxlevels]]$result[taxa, "MRCA", drop=FALSE]  # extract MRCA column out of result table for specified taxa only
          allnodes[[namenod]] <- tmp  # add the extracted table as sub-entry to output list, named according to taxlevel
        }
    }
    if (length(allnodes) == 0) {
      stop('No results found. Consider reviewing taxlevels argument.')
    }
    allnodes  # return final list with MRCA nodes
}
