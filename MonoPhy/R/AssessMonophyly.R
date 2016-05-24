# assesses monophyly of genera (or customized units) and makes the result available in different ways (tables, onjects, plot...)
# written by Orlando Schwery 2015
AssessMonophyly <-
function(tree, taxonomy=NULL, verbosity=5, outliercheck=TRUE, outlierlevel=0.5, taxizelevel= NULL, taxizedb='ncbi', taxizepref='ncbi', taxask=FALSE, taxverbose=FALSE) {
# initial tests and data preparation
    if (!is.rooted(tree)) {  # checks and returns error if tree is not rooted
        stop('Phylogeny must be rooted!')
    }
    if (is.null(taxonomy)) {  # if no taxonomy file is specified, extract list of genera from tree's tip labels
        for (i in 1:length(tree$tip.label)) {  # loop through tip labels of tree
            if (grepl(("_| "), tree$tip.label[i]) == FALSE) {  # checks if genus and species epithet of tip labels are separated by space or underscore and returns error if not
                stop('Tip labels do not contain underscore separating genus name from species epithet!')
            }
        }
        f <- function(s) strsplit(s, ("_| "))[[1]][1]  # function with split criteria: split names at underscore and keep first part (genus name)
        split.taxa <- sapply(tree$tip.label, f)  # apply split functon to tree
        taxa <- as.vector(unique(split.taxa))  # create vector of genera in tree without duplicates
        taxsets <- c('taxa')  # assign 'taxa' as taxsetnames
    } else {
        if (!is.null(taxonomy) && taxonomy != 'taxize') {  # if argument 'taxonomy' is not NULL and not taxize, use loaded taxonomy file
            if (length(taxonomy[, 1]) != length(tree$tip.label)) {  # checks and returns error if taxonomy file has more or less entries than tree has tips
                stop('Number of rows of taxonomy file is not equal to number of taxa (note: if your table has a header, you should specify header=TRUE when importing)!')
            }
            if (length(taxonomy[1, ]) < 2) {  # checks and returns error if taxonomy file doesn't have at least two columns
                stop('Taxonomy file needs at least 2 columns: tip labels and taxonomic group!')
            }
            taxchecktree <- c()  # create empty vector to be filled with presence of tip labels in taxfile
            for (itaxcheck in 1:length(tree$tip.label)) {  # loop through tip labels
                taxchecktree <- c(taxchecktree, tree$tip.label[itaxcheck] %in% taxonomy[, 1])  # check for every name in tree if it is present in taxonomy file
            }
            taxintruderstree <- c()  # create empty vector to be filled with tip labels that are abesent in taxfile
            if ('FALSE' %in% taxchecktree) {  # if there are any missing names in taxonomy file...
                positionstree <- grep('FALSE', taxchecktree)  # ...get their position...
                taxintruderstree <- c(taxintruderstree, tree$tip.label[positionstree])  # ...add tip labels of that position to vector
                message(paste('\n'), appendLF=TRUE)
                message(paste('Tip-labels which do not occur in taxonomy file:', '[', length(taxintruderstree), '/', length(tree$tip.label), ']', collapse=" "), appendLF=TRUE)  # display numbers of missing tip labels in taxonomy file
                message(paste(taxintruderstree, collapse=", "), appendLF=TRUE)  # display names of missing tips
                message(paste('\n'), appendLF=TRUE)
            }
            taxcheckfile <- c()  # create empty vector to be filled with presence of taxfile names in tip labels
            for (itaxcheck2 in 1:length(taxonomy[, 1])) {  # loop through names column of taxfile
                taxcheckfile <- c(taxcheckfile, taxonomy[itaxcheck2, 1] %in% tree$tip.label)  # check for every name in taxonomy file if it is present in tip labels of the tree
            }
            taxintrudersfile <- c()  # create empty vector to be filled with taxfile names that are abesent in tip labels
            if ('FALSE' %in% taxcheckfile) {  # if there are any missing names in tip labels...
                positionsfile <- grep('FALSE', taxcheckfile)  # ...get their position...
                taxintrudersfile <- c(taxintrudersfile, as.character(taxonomy[positionsfile, 1]))  # ...add names of that position to vector
                message(paste('Taxon names in file which do not occur in tip-labels of tree:', '[', length(taxintrudersfile), '/', length(taxonomy[, 1]), ']', collapse=" "), appendLF=TRUE)  # display numbers of missing names in tip labels
                message(paste(taxintrudersfile, collapse=", "), appendLF=TRUE)  # display names of missing taxa
                message(paste('\n'), appendLF=TRUE)
            }
            if ('FALSE' %in% (taxchecktree)) {  # if missing names in file, stop and display error
                stop('The taxon names of tree and taxonfile do not match (see above)!')
            }
            if ('FALSE' %in% (taxcheckfile)) {  # if missing names in tree, stop and display error
                stop('The taxon names of tree and taxonfile do not match (see above)!')
            }
            taxsets <- list()  # create empty list to fill with taxonomic units
            for (jtax in 1:(length(taxonomy[1, ]) - 1)) {  # loop through taxon file
                nametax <- paste("taxa", jtax, sep="")  # create taxon name label
                tmp <- as.vector(unique(taxonomy[, (jtax) + 1]))  # if all is correct, makes vector of taxonomic units (without doubles)
                taxsets[[nametax]] <- tmp  # add names vector to taxets list, labelled with taxon name label
            }
        } else {  # if not NULL and not taxfile but taxize
            if (taxonomy == 'taxize') {  # build taxonomy file from web ressources using taxize
                taxafromweb <- tax_name(tree$tip.label, get=taxizelevel, db=taxizedb, pref=taxizepref, ask=taxask, verbose=taxverbose)  # get taxonomy data from web
                taxafromwebtable <- matrix(data = NA, nrow=length(tree$tip.label),ncol=ncol(taxafromweb) - 1)  #build empty matrix with dimensions by number of tips and retrieved taxonomic data
                taxafromwebtable[, 1] <- tree$tip.label  # add tip names from tree
                if (length(unique(taxafromweb[, 3])) == 1 & is.na(unique(taxafromweb[, 3])[1])) {  # check if there was any information retrieved and display error if not
                    stop('There was no data found for any of the tips!')
                }
                if (nrow(taxafromweb) > length(tree$tip.label)) {  # check if more record entires retrieved than tips in tree
                    taxafromweb <- unique(taxafromweb[, 2:3])  # if entires from two databases are the same, delete one
                    rownames(taxafromweb) <- c(1:nrow(taxafromweb))  # renumber rows
                }
                if (nrow(taxafromweb) > length(tree$tip.label)) {  # check again if more record entires retrieved than tips in tree
                    temp <- taxafromweb  # copy table to temporary object
                    droppers <- c()  # create empty vector to be filled with entires to be dropped
                    counterweb <- table(temp[, 1], useNA='ifany')  # create counting table with number of entries per taxon in retrieved data table
                    for (iwebtax in 1:nrow(temp)) {  # loop through retrieved table
                        if (counterweb[temp[iwebtax, 1]] > 1) {  # check if more than one entry per species is retrieved (according to count table)
                            if (is.na(temp[iwebtax, 2])) {  # if the duplicate currently being looked at did not retrieve a result...
                                droppers <- c(droppers, iwebtax)  # ...add it to droplist
                            }
                        }
                        temp2 <- temp[-droppers, ]  # remove drop entires
                        rownames(temp2) <- c(1:nrow(temp2))  # renumber rows
                        taxafromweb <- temp2  # replace retrieved table with cleaned up version
                    }
                }
                if (colnames(taxafromweb[1]) == "db") {  # if the first column is the name of the database (which is getting lost if duplicates are removed)...
                    taxafromweb <- taxafromweb[, -1]  # ...that first column will be dropped to get the same table format.
                }
                for (iweb in 2:(ncol(taxafromweb))) {  #add acquired taxon names for tips for each taxonomic level acquired
                    taxafromwebtable[, iweb] <- taxafromweb[, iweb]  # add retrieved entries to matrix
                    rownames(taxafromweb) <- c(1:nrow(taxafromweb))  # renumber row names
                }
                taxafromwebtable[is.na(taxafromwebtable)] <- "unknown"  # replace all NAs with "unknown"
                taxonomy <- as.data.frame(taxafromwebtable)  # turn matrix into data frame and feed to further function
                taxsets <- list()  # create empty list to be filled with taxonomy lists
                for (jtax in 1:(length(taxonomy[1, ]) - 1)) {  # loop through taxonomy file
                    nametax <- paste("taxa", jtax, sep="")  # create taxon name label
                    tmp <- as.vector(unique(taxonomy[, (jtax) + 1]))  # if all is correct, makes vector taxonomic units (without doubles)
                    taxsets[[nametax]] <- tmp  # add names vector to taxets list, labelled with taxon name label
                }
            }
        }
    }
# actual assessment
    finallist <- list()  # create empty list to be filled with final results
    for (ifullround in 1:length(taxsets)) {  # Assess monophyly for every taxon set used
        if (is.null(taxonomy)) {  # if no taxonomy file specified...
            taxa <- taxa  # ...assign taxon list
        } else {  # if taxonomy file or taxize...
            taxa <- unlist(taxsets[ifullround])  # ...assign taxon list of current taxlevel...
            taxa <- taxa[!taxa %in% c("unknown", NA)]  # ...and remove NA's and unknown taxa
        }
    # create empty objects to be filled by function
        intruder.genus <- list()  # empty list for genera causing non-monophyly
        intruder.genus.full <- c()  # empty vector for ALL general causing non-monophyly
        intruder.species <- list()  # empty list for species causing non-monophyly
        intruder.species.full <- c()  # empty vector for ALL species causing non-monophyly
        intruder.names <- c()  # empty vector for names for intruder sub-lists
        if (outliercheck == TRUE) {  # if outliers should be checked for (add additional objects and columns/rows
            outlist.summary <- matrix(NA, nrow=6, ncol=3)  # create final output summary matrix
            dfheaders <- c("Taxon", "Monophyly", "MRCA", "#Tips", "Delta-Tips", "#Intruders", "Intruders", "#Outliers", "Outliers")  # headers for 'outlist'
            outlist <- matrix(NA, nrow=length(taxa), ncol=9)  # final output matrix
            outlier.species <- list()  # list of species causing non-monophyly as outliers
            outlier.species.full <- c()  # vector of ALL species causing non-monophyly as outliers
            outlier.names <- c()  # names for outlier sub-lists
        } else {  # if no outliers being checked for
            outlist.summary <- matrix(NA, nrow=5, ncol=3)  # final output summary matrix
            dfheaders <- c("Taxon", "Monophyly", "MRCA", "#Tips", "Delta-Tips", "#Intruders", "Intruders")  # headers for 'outlist'
            outlist <- matrix(NA, nrow=length(taxa), ncol=7)  # final output matrix
        }
        tip.states.matrix <- matrix(NA, nrow=length(tree$tip.label), ncol=3)  # states for plotting
# loop assessing monophyly
        for (i in 1:length(taxa)) {  # loop through every genus in the tree
            if (is.null(taxonomy)) {  # genera extracted from tip labels if no taxonomy file loaded
                ancnode <- getMRCA(tree, tip=c(tree$tip.label[c(grep(paste("^", taxa[i], "_", sep=""), tree$tip.label))]))  # determine Most Recent Common Ancestor for all taxa of genus
            } else {  # units taken from file if loaded
                subtips <- subset(taxonomy, as.character(taxonomy[, (ifullround+1)]) == as.character(taxa[i]))  # get tips associated with current group
                subtipsnr <- c()  # create empty vector to be filled with tip numbers
                for (sbts in 1: nrow(subtips)) {  # looping through all subtips assigned to current group
                    sbtname <- subtips[sbts, 1]  # extract name
                    sbtnr <- which(tree$tip.label == sbtname)  # extract number of subtip
                    subtipsnr <- c(subtipsnr, sbtnr)  # add subtip nr to numbers vector
                }
                ancnode <- getMRCA(tree, tip=c(subtipsnr))  # get MRCA for tips associated with group
            }
            if (length(ancnode) == 0) { # if monotypic i.e. only tip of given group
                if (outliercheck == TRUE) {  # if outliers are checked for (more columns)
                    outlist[i, ] <- c(taxa[i], "Monotypic", "NA", 1, "NA", "NA", "", "NA", "")  # UPDATE OUTPUT MATRIX, mark as monotypic if only one tip for this genus
                } else {  # if outliers are not checked for (less columns)
                    outlist[i, ] <- c(taxa[i], "Monotypic", "NA", 1, "NA", "NA", "")  # UPDATE OUTPUT MATRIX, mark as monotypic if only one tip for this genus
                }
            } else {  # if not monotypic
                anctips <- getDescendants(tree, ancnode)  # determine all descendants of previously determined MRCA
                ancnames <- tree$tip.label[c(anctips)]  # extract names of those descendants
                ancnames <- ancnames[!is.na(ancnames)]  # ommit NA's (caused by descendants which are internal nodes and not tips)
                if (is.null(taxonomy)) {  # genera extracted from tip labels if no taxonomy file loaded
                    taxtips <- tree$tip.label[c(grep(paste("^", taxa[i], "_", sep=""), tree$tip.label))]  # get tip names of genus in question
                } else {  # if taxonomy file loaded
                    taxtips <- subtips[, 1]  # get vector of tip names of genus in question
                }
                if (length(ancnames) == length(taxtips)) {  # determine if all MRCA descendants = genus members. Genus is monophyletic if yes.
                    if (outliercheck == TRUE) {  # if outliers are checked for (more columns)
                        outlist[i, ] <- c(taxa[i], "Yes", ancnode, length(taxtips), "0", "0", "", "NA", "")  # UPDATE OUTPUT MATRIX, mark as monophyletic
                    } else {  # if outliers are checked for (less columns)
                        outlist[i, ] <- c(taxa[i], "Yes", ancnode, length(taxtips), "0", "0", "")  # UPDATE OUTPUT MATRIX, mark as monophyletic
                    }
                } else {  # if taxon is not monophyletic
                    intruder.tips <- setdiff(ancnames, taxtips)  # determine intruders tip labels, i.e. descendants of MRCA which are not genus members
                    if (is.null(taxonomy)) {  # get intruder genus names if no taxonomy file loaded
                        f2 <- function(s) strsplit(s, ("_| "))[[1]][1]  # function with split criteria: split names at underscore and keep first part (genus name)
                        split.taxa2 <- sapply(intruder.tips, f2)  # apply split function to intruder tip labels
                        intruder.taxa <- as.vector(unique(split.taxa2))  # create vector of intruder genera
                    } else {  # get intruder taxonomic levels from taxonomy file
                        subtaxa <- c()  # create empty vector to be filled with intruder taxnames
                        for (j in 1:length(intruder.tips)) {  # loop through intruder tips
                            subtaxon <- rbind(subset(taxonomy, taxonomy[, 1] == intruder.tips[j]))  # extract taxon for each intruder tip...
                            subtaxa <- rbind(subtaxa, subtaxon)  # ... and add them up
                        }
                        intruder.taxa <- as.vector(unique(subtaxa[, ifullround + 1]))  # create vector of intruder taxa
                    }
                    outlier.tips <- c()  # create empty vector to be filled with outlier tips
                    if (outliercheck == TRUE) {  # distinguish outliers if TRUE
                      if (length(Children(tree, ancnode)) > 2) {
                        outlier.tips <- c()
                      } else {
                        tiplevels <- length(taxtips) / length(ancnames)  # determine fraction of total tips to members of focal taxon among descendants of current MRCA
                        if (tiplevels < outlierlevel ) {  # check if meeting criteria
                            start.node <- ancnode  # set MRCA node as starting point
                            while (tiplevels < outlierlevel) {  # search for subclade that meets criteria
                                subtaxtips <- c()  # reset taxtips
                                subancnames <- c()  # reset ancnames
                                parent.node <- start.node # set parent node
                                daughter.nodes <- Children(tree, parent.node) # find direct descendant nodes
                                if (length(daughter.nodes) > 2) {  # if multifurcation, quit
                                  anctips1 <- getDescendants(tree, parent.node)
                                  ancnames1 <- tree$tip.label[c(anctips1)]  # extract names of those descendants
                                  subancnames <- ancnames1[!is.na(ancnames1)]  # ommit NA's (caused by descendants which are internal nodes and not tips)
                                  subtaxtips <- intersect(taxtips, subancnames)
                                  start.node <- parent.node
                                    break # end loop
                                }
                                daughter1 <- daughter.nodes[1]  # assign daughter node 1 separately
                                daughter2 <- daughter.nodes[2]  # assign daughter node 2 separately
                            # prepare descendants of daughter 1
                                anctips1 <- getDescendants(tree, daughter1)  # determine all descendants of daughter1
                                ancnames1 <- tree$tip.label[c(anctips1)]  # extract names of those descendants
                                ancnames1 <- ancnames1[!is.na(ancnames1)]  # ommit NA's (caused by descendants which are internal nodes and not tips)
                                taxtips1 <- intersect(taxtips, ancnames1)  # get taxon members of subclade1
                            # prepare descendants of daughter2
                                anctips2 <- getDescendants(tree, daughter2)  # determine all descendants of daughter2
                                ancnames2 <- tree$tip.label[c(anctips2)]  # extract names of those descendants
                                ancnames2 <- ancnames2[!is.na(ancnames2)]  # ommit NA's (caused by descendants which are internal nodes and not tips)
                                taxtips2 <- intersect(taxtips, ancnames2)  # get taxon members of subclade2
                            # pick from the daughter nodes
                                nodechoice <- which(c(length(taxtips1), length(taxtips2)) == max(c(length(taxtips1), length(taxtips2))))  # determine daughter with more tips of focal taxon
                                if (length(nodechoice) > 1) {  # if equal number of taxon tips in both daughers:
                                    nodechoice2 <- which(c((length(taxtips1) / length(anctips1)), (length(taxtips2) / length(anctips2))) == max(c((length(taxtips1) / length(anctips1)), (length(taxtips2) / length(anctips2)))))  # determine daugther with higher ratio of focal taxon
                                    if (length(nodechoice2) > 1) {  # if equal ratio of taxon tips in both daughters: keep both daughters as core clade
                                        subtaxtips <- c(taxtips1, taxtips2)  # merge their taxon members
                                        subancnames <- c(ancnames1, ancnames2)  # merge their node descendants
                                        start.node <- parent.node  # set their parent node as start point
                                          break  # end loop
                                    } else if (nodechoice2 == 1) {  # if daughter 1 chosen: set as new start-node and -clade
                                        subtaxtips <- taxtips1  # assign taxon members
                                        subancnames <- ancnames1  # assign node descendants
                                        start.node <- daughter1  # set node as start point
                                    } else if (nodechoice2 == 2) {  # if daughter 2 chosen: set as new start-node and -clade
                                        subtaxtips <- taxtips2  # assign taxon members
                                        subancnames <- ancnames2  # assign node descendants
                                        start.node <- daughter2  # set node as start point
                                    }
                                } else if (nodechoice == 1) {  # if daughter 1 chosen: set as new start-node and -clade
                                    subtaxtips <- taxtips1  # assign taxon members
                                    subancnames <- ancnames1  # assign node descendants
                                    start.node <- daughter1  # set node as start point
                                } else if (nodechoice == 2) {  # if daughter 2 chosen: set as new start-node and -clade
                                    subtaxtips <- taxtips2  # assign taxon members
                                    subancnames <- ancnames2  # assign node descendants
                                    start.node <- daughter2  # set node as start point
                                }
                                tiplevels <- length(subtaxtips) / length(subancnames)  # reassess status of current clade
                            }
                            if (tiplevels < 1 & length(daughter.nodes) <=2) {  # if intruders are present after outliercheck, check if early-diverging
                                EDtaxtips1 <- c()  # create empty vector to be filled with early diverging tips
                                EDtaxtips2 <- c()  # create empty vector to be filled with early diverging tips
                            # search for node whose daughers both include members of the focal taxon
                                repeat {  # continue search until break conditions are met
                                    EDparent.node <- start.node # set parent node
                                    if (EDparent.node <= length(tree$tip.label)) {  # if parent node is actually a tip
                                        subtaxtips <- tree$tip.label[EDparent.node]  # assign that tip to taxon members
                                        subancnames <-tree$tip.label[EDparent.node]  # assign that tip to node descendants
                                        start.node <- EDparent.node  # set that node as start point
                                        break  # end checking for early diverging
                                    }
                                    EDdaughter.nodes <- Children(tree, EDparent.node)  # find direct descendant nodes
                                    EDdaughter1 <- EDdaughter.nodes[1]  # assign daughter node 1 separately
                                    EDdaughter2 <- EDdaughter.nodes[2]  # assign daughter node 2 separately
                                # prepare descendants of daughter 1
                                    EDanctips1 <- getDescendants(tree, EDdaughter1)  # determine all descendants of EDdaughter1
                                    EDancnames1 <- tree$tip.label[c(EDanctips1)]  # extract names of those descendants
                                    EDancnames1 <- EDancnames1[!is.na(EDancnames1)]  # ommit NA's (caused by descendants which are internal nodes and not tips)
                                    EDtaxtips1 <- intersect(taxtips, EDancnames1)  # get taxon members of subclade1
                                # prepare descendants of daughter 2
                                    EDanctips2 <- getDescendants(tree, EDdaughter2)  # determine all descendants of EDdaughter2
                                    EDancnames2 <- tree$tip.label[c(EDanctips2)]  # extract names of those descendants
                                    EDancnames2 <- EDancnames2[!is.na(EDancnames2)]  # ommit NA's (caused by descendants which are internal nodes and not tips)
                                    EDtaxtips2 <- intersect(taxtips, EDancnames2)  # get taxon members of subclade2
                                # select node to continue
                                    if (length(EDtaxtips1) == 0) {  # if EDdaugher1 has no decendant of focal taxon, continue with EDdaughter2
                                        start.node <- EDdaughter2  # set daughter 2 as new start node
                                    } else if (length(EDtaxtips2) == 0) {  # if EDdaugher2 has no decendant of focal taxon, continue with EDdaughter1
                                        start.node <- EDdaughter1  # set daughter 1 as new start node
                                    } else {  # if both daughers have descendants of focal taxon, stop by keeping their parent node as node of core clade
                                        subtaxtips <- c(EDtaxtips1, EDtaxtips2)  # merge their taxon members
                                        subancnames <- c(EDancnames1, EDancnames2)  # merge their node descendants
                                        start.node <- EDparent.node # set their parent node as start point
                                        break  # end checking for early diverging
                                    }
                                }
                            }
                            outlier.tips <- setdiff(taxtips, subtaxtips)  # determine outliers
                            if (length(outlier.tips) != 0) {  # if there ARE outliers...
                                outlier.species <- c(outlier.species, list(Tips=outlier.tips))  # ...update list of outlier tip labels
                                outlier.species.full <- c(outlier.species.full, outlier.tips)  # ...update vector of ALL outlier species
                                outlier.names <- c(outlier.names, taxa[i])  # ...update names vector for outliers
                            }
                            intruder.tips <- setdiff(subancnames, subtaxtips)  # determine intruders
                            if (is.null(taxonomy)) {  # get intruder genus names if no taxonomy file loaded
                                f2 <- function(s) strsplit(s, ("_| "))[[1]][1]  # function with split criteria: split names at underscore and keep first part (genus name)
                                split.taxa2 <- sapply(intruder.tips, f2)  # apply split function to intruder tip labels
                                intruder.taxa <- as.vector(unique(split.taxa2))  # create vector of intruder genera
                            } else {  # get intruder taxonomic levels from taxonomy file
                                subtaxa <- c()  # create empty vector to be filled with taxon names
                                for (j in 1:length(intruder.tips)) {  # loop through intruder tips
                                    subtaxon <- rbind(subset(taxonomy, taxonomy[, 1] == intruder.tips[j]))  # extract taxon for each intruder tip...
                                    subtaxa <- rbind(subtaxa, subtaxon)  # ... and add them up
                                }
                                intruder.taxa <- as.vector(unique(subtaxa[, ifullround + 1]))  # create vector of intruder taxa (without doubles)
                            }
                        }
                    }
                    if (length(intruder.taxa) != 0) {  # if there ARE intruders...
                        intruder.genus <- c(intruder.genus, list(Taxa=intruder.taxa))  # ...update list of intruder genera
                        intruder.genus.full <- c(intruder.genus.full, intruder.taxa)  # ...update vector of ALL intruder genera
                        intruder.species <- c(intruder.species, list(Tips=intruder.tips))  # ...update list of intruder tip labels
                        intruder.species.full <- c(intruder.species.full, intruder.tips)  # ...update vector of ALL intruder species
                        intruder.names <- c(intruder.names, taxa[i])  # ...update names vector for intruders
                    }
                    if (outliercheck == TRUE) {  # if outliers are being checked for
                        if (length(intruder.taxa) <= verbosity) {  # if less intruding genera than specified in verbosity...
                            if (length(outlier.tips) <= verbosity) {  # if less intruding tips than specified in verbosity...
                                outlist[i, ] <- c(taxa[i], "No", ancnode, length(taxtips), (length(ancnames) - length(taxtips)), length(intruder.taxa), paste(intruder.taxa, collapse=", "), length(outlier.tips), paste(outlier.tips, collapse=", "))  # UPDATE OUTPUT MATRIX, mark as non-monophyletic and list intruder genera
                            } else {  # if more intruder tips than verbosity level...
                               outlist[i, ] <- c(taxa[i], "No", ancnode, length(taxtips), (length(ancnames) - length(taxtips)), length(intruder.taxa), paste(intruder.taxa, collapse=", "), length(outlier.tips), paste(outlier.tips[1], "and", (length(outlier.tips) - 1), "more.", collapse=", "))  # UPDATE OUTPUT MATRIX, mark as non-monophyletic and list intruder genera
                            }
                        } else {  # if more intruder taxa than verbosity level...
                            if (length(outlier.tips) <= verbosity) {  # if less intruding tips than specified in verbosity...
                                outlist[i, ] <- c(taxa[i], "No", ancnode, length(taxtips), (length(ancnames) - length(taxtips)), length(intruder.taxa), paste(intruder.taxa[1], "and", (length(intruder.taxa) - 1), "more.", collapse=", "), length(outlier.tips), paste(outlier.tips, collapse=", "))  # UPDATE OUTPUT MATRIX, mark as non-monophyletic and list intruder genera
                            } else {  # if more intruder tips than verbosity level...
                                outlist[i, ] <- c(taxa[i], "No", ancnode, length(taxtips), (length(ancnames) - length(taxtips)), length(intruder.taxa), paste(intruder.taxa[1], "and", (length(intruder.taxa) - 1), "more.", collapse=", "), length(outlier.tips), paste(outlier.tips[1], "and", (length(outlier.tips) - 1), "more.", collapse=", "))  # UPDATE OUTPUT MATRIX, mark as non-monophyletic and list intruder genera
                            }
                        }
                    } else {  # if outliers are not checked for...
                        if (length(intruder.taxa) <= verbosity) {  # if less intruding genera than specified in verbosity...
                            outlist[i, ] <- c(taxa[i], "No", ancnode, length(taxtips), (length(ancnames) - length(taxtips)), length(intruder.taxa), paste(intruder.taxa, collapse=", "))  # UPDATE OUTPUT MATRIX, mark as non-monophyletic and list intruder genera
                        } else {  # if more intruder taxa than verbosity level...
                            outlist[i, ] <- c(taxa[i], "No", ancnode, length(taxtips), (length(ancnames) - length(taxtips)), length(intruder.taxa), paste(intruder.taxa[1], "and", (length(intruder.taxa) - 1), "more.", collapse=", "))  # UPDATE OUTPUT MATRIX, mark as non-monophyletic and list intruder genera
                        }
                    }
                  }
                }
            }
        }
# prepare outputs
        intruder.genus.all <- unique(intruder.genus.full)  # vector of ALL intruder genera without doubles
        intruder.species.all <- unique(intruder.species.full)  # vector of ALL intruder species without doubles
        outlier.species.all <- c()  # create empty vector to be filled with outlier species names
        if (outliercheck == TRUE) {  # if outliers are being checked for...
            outlier.species.all <- unique(outlier.species.full)  # ...create vector of ALL outlier species without doubles
        }
        outframe <- data.frame(outlist)  # turn final output matrix into data frame
        names(outframe) <- dfheaders  # apply headers
        rownames(outframe) <- outframe[, 1]  # assign first column as row names
        outframe[, 1] <- NULL  # delete first column (since now row names)
        names(intruder.genus) <- intruder.names  # apply names (of intruded genera) to intruder genus lists
        names(intruder.species) <- intruder.names  # apply names (of intruded genera) to intruder species lists
        if (is.null(taxonomy)) {  # if no taxonomy is specified (i.e. genus names are grepped from tip labels
            tip.states.matrix[, 1] <- tree$tip.label  # adding species names to matrix
            f3 <- function(s) strsplit(s, ("_| "))[[1]][1]  # function with split criteria: split names at underscore and keep first part (genus name)
            tip.states.matrix[, 2] <- sapply(tip.states.matrix[, 1], f3)  # apply split function, create genus name column
        } else {  # if taxonomy file was specified (or created using taxize)
            tip.states.matrix[, 1] <- as.vector(taxonomy[, 1])  # adding species names to matrix
            tip.states.matrix[, 2] <- as.vector(taxonomy[, ifullround + 1])  # create taxon name column
        }
        for (i in 1:length(tree$tip.label)) {  # loop through tip labels
            if (tip.states.matrix[i, 1] %in% intruder.species.all == TRUE) {  # ...if species is in global intruder list...
                tip.states.matrix[i, 3] <- "Intruder"  # ...score species as intruder
            } else if (tip.states.matrix[i, 1] %in% outlier.species.all == TRUE) {  # ...if species is in global outlier list...
                tip.states.matrix[i, 3] <- "Outlier"  # ...score species as outlier
            } else if (outframe[tip.states.matrix[i, 2], "Monophyly"] == "Monotypic") {  # ...if  species is monotypic...
                tip.states.matrix[i, 3] <- "Monophyletic"  # ... score as monophyletic
            } else if (outframe[tip.states.matrix[i, 2], "Monophyly"] == "Yes"){  # ... if species is monophyletic...
                tip.states.matrix[i, 3] <- "Monophyletic"  # ... score as monophyletic
            } else if (outframe[tip.states.matrix[i, 2], "Monophyly"] == "No") {  # ... if species is not monophyletic...
                tip.states.matrix[i, 3] <- "Non-Monophyletic"  # ... score as non-monophyletic
            } else if (tip.states.matrix[i, 2] == "unknown" | is.na(tip.states.matrix[i, 2])) {  # if assignment to taxonomic group is unknown or NA...
                tip.states.matrix[i, 3] <- "unknown"  # ...set monophyly state to 'unknown'
            }
        }
        tip.states.frame <- as.data.frame(tip.states.matrix)  # turn tip states matrix into data frame
        colnames(tip.states.frame) <- c("Tip", "Taxon", "Status")  # assign column names
        if (outliercheck == TRUE) {  # if outliers are being checked for...
            names(outlier.species) <- outlier.names  # ...apply names (of outlier genera) to outlier species lists
            outlist.summary[, 1] <- c("Total", "Monophyletic", "Non-Monophyletic", "Monotypic", "Intruder", "Outlier")  # ...add row names for summary table
        } else {  # if outliers are not checked for...
            outlist.summary[, 1] <- c("Total", "Monophyletic", "Non-Monophyletic", "Monotypic", "Intruder")  # ...add row names for summary table
        }
    # make counts for summary table
        counttable <- table(outframe[, "Monophyly"])  # tabulate monophyly results
        countframe <- as.data.frame(counttable)  # turn into data frame
        rownames(countframe) <- countframe[, 1]  # assign first column as row names
        countframe[, 1] <- NULL  # delete first column (since now row names)
    # make tip counts
        mono.count <- c()  # create empty vector to be filled with numbers of monophyletic tips
        nonmono.count <- c()  # create empty vector to be filled with numbers of non-monophyletic tips
        for (itaxcount1 in 1:length(taxa)) {  # for every taxon...
            if (outlist[itaxcount1, 2] == "Yes") {  # ...if monophyly is inferred to be true...
                mono.count <- c(mono.count, as.numeric(outlist[itaxcount1, 4]))  # ...add the number of tips for that taxon to the monophyly counter
            }
        }
        for (itaxcount2 in 1:length(taxa)) {  # for every taxon...
            if (outlist[itaxcount2,2] == "No") {  # ...if monophyly is inferred to be not true...
                nonmono.count <- c(nonmono.count, as.numeric(outlist[itaxcount2, 4]))  # ...add the number of tips for that taxon to the non-monophyly counter
            }
        }
        if (outliercheck == TRUE) {  # if outliers are being checked for...
            outlist.summary[, 2] <- c(length(taxa), countframe["Yes", "Freq"], countframe["No", "Freq"], countframe["Monotypic", "Freq"], length(intruder.genus.all), length(outlier.names))  # populate summary table with counts of respective groups (taxa)
            outlist.summary[, 3] <- c(length(tree$tip.label), sum(mono.count), sum(nonmono.count), countframe["Monotypic", "Freq"], length(intruder.species.all), length(outlier.species.all))  # populate summary table with counts of respective groups (tips)
        } else {  # ...if outliers are not being checked for...
            outlist.summary[, 2] <- c(length(taxa), countframe["Yes", "Freq"], countframe["No", "Freq"], countframe["Monotypic", "Freq"], length(intruder.genus.all))  # populate summary table with counts of respective groups (taxa)
            outlist.summary[, 3] <- c(length(tree$tip.label), sum(mono.count), sum(nonmono.count), countframe["Monotypic", "Freq"], length(intruder.species.all))  # populate summary table with counts of respective groups (tips)
        }
        outframe.summary <- data.frame(outlist.summary)  # turn final output matrix summary into data frame
        rownames(outframe.summary) <- outframe.summary[, 1]  # assign first column as row names
        outframe.summary[, 1] <- NULL  # delete first column (since now row names)
        colnames(outframe.summary) <- c("Taxa", "Tips")  # assign column names to summary data frame
        if (outliercheck == TRUE) {  # if outliers are being checked for...
            outputlist <- list(IntruderTaxa=intruder.genus, IntruderTips=intruder.species, OutlierTaxa=outlier.names, OutlierTips=outlier.species, result=outframe, summary=outframe.summary, TipStates=tip.states.frame) # concatenate intruder lists, final output table and summmary to one list object
        } else {  # ...if outliers are not being checked for...
            outputlist <- list(IntruderTaxa=intruder.genus, IntruderTips=intruder.species, result=outframe, summary=outframe.summary, TipStates=tip.states.frame) # concatenate intruder lists, final output table and summmary to one list object
        }
        if (is.null(taxonomy)) {  #Add taxlevel names if no input taxonomy...
          nameout <- "Genera"  # ...name level "Genera"
        } else if (!is.null(taxonomy) && is.null(taxizelevel)) {  # if taxtable used...
          if (colnames(taxonomy)[ifullround] == paste("V", ifullround, sep="")) {  # if taxonomy table doesn't have header (or it is empty)...
            nameout <- paste("Taxlevel", ifullround, sep="_")  # ...name level of output subsection 'taxlevel#'
          } else {  # if taxonomy table has header...
            nameout <- colnames(taxonomy)[ifullround+1]  # ...name level accroding to header entry
          }
        } else if (!is.null(taxizelevel)) {  # if taxonomy from taxize...
            nameout <- taxizelevel[ifullround]  # ...name level according to taxize query
        }
        finallist[[nameout]] <- outputlist #add list for this round of the loop to final list
    }
    finallist  # return final outputlist
}
