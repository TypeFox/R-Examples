

###############################
#### SUBSETTING PROCEDURES ####
###############################

#########################
## obkSequences method ##
#########################
setMethod("subset", "obkSequences", function(x, sequences=NULL, locus=NULL, individuals=NULL,
                                        date.from=NULL, date.to=NULL, date.format=NULL, ...){
    ## ESCAPE IF NOTHING IN THE DATASET ##
    if(get.nsequences(x)==0) return(x)

    ## SUBSET BY SEQUENCE ID ##
    if(!is.null(sequences)){
        ## handle non-character argument ##
        if(!is.character(sequences)){
            sequences <- get.sequences(x)[sequences]
        }

        ## check that all sequences are known ##
        if(!all(sequences %in% get.sequences(x))){
            sequences <- na.omit(sequences)
            temp <- paste(sequences[!sequences %in% get.sequences(x)], collapse=", ")
            warning(paste("The following sequences were not found in the data:", temp))
            sequences <- sequences[sequences %in% get.sequences(x)]
        }

        ## subset x@dna ##
        x@dna <- get.dna(x, id=sequences)

        ## subset x@meta ##
        x@meta <- x@meta[sequences,,drop=FALSE]

    } # end subset by sequence ID


    ## SUBSET BY LOCUS ##
    if(!is.null(locus)){
        locus <- na.omit(locus)

        ## handle non-character argument ##
        if(!is.character(locus)){
            locus <- get.locus(x)[locus]
        }

        ## check that all loci are known ##
        if(!all(locus %in% get.locus(x))){
            temp <- paste(locus[!locus %in% get.locus(x)], collapse=", ")
            warning(paste("The following loci were not found in the data:", temp))
            locus <- locus[locus %in% get.locus(x)]
        }

        ## subset x@dna ##
        x@dna <- x@dna[locus]
        seq.tokeep <- unlist(lapply(x@dna, rownames))

        ## subset x@meta ##
        x@meta <- x@meta[seq.tokeep,,drop=FALSE]

    } # end subset by locus


     ## SUBSET BY INDIVIDUAL ID ##
    if(!is.null(individuals)){
        individuals <- na.omit(individuals)
        ## handle non-character argument ##
        if(!is.character(individuals)){
            individuals <- get.individuals(x)[individuals]
        }

        ## check that all individuals are known ##
        if(!all(individuals %in% get.individuals(x))){
            temp <- paste(individuals[!individuals %in% get.individuals(x)], collapse=", ")
            warning(paste("The following individuals were not found in the data:", temp))
            individuals <- individuals[individuals %in% get.individuals(x)]
        }

        ## find sequences to retain ##
        seq.tokeep <- rownames(x@meta)[x@meta$individualID %in% individuals]

        ## subset by sequences ##
        x <- subset(x, sequences=seq.tokeep)
    } # end subset by individual ID


    ## DATES FROM ... ##
    if(!is.null(date.from)){
        date.from <- na.omit(date.from)

        ## process date ##
        date.from <- .process.Date(date.from, format=date.format)

        ## find sequences to retain ##
        seq.tokeep <- x@meta$date >= date.from

        ## subset by sequences ##
        x <- subset(x, sequences=seq.tokeep)
    } # end subset by date.from


    ## DATES TO ... ##
    if(!is.null(date.to)){
        date.to <- na.omit(date.to)

        ## process date ##
        date.to <- .process.Date(date.to, format=date.format)

        ## find sequences to retain ##
        seq.tokeep <- x@meta$date <= date.to

        ## subset by sequences ##
        x <- subset(x, sequences=seq.tokeep)
    } # end subset by date.to

    return(x)
}) # end subset for obkSequences






########################
## obkContacts method ##
########################
setMethod("subset", "obkContacts", function(x, individuals=NULL, date.from=NULL, date.to=NULL,
                                            date.format=NULL, ...){
    ## SUBSET BY INDIVIDUALS ##
    if(!is.null(individuals)){
        individuals <- na.omit(individuals)

        ## handle non-character argument ##
        if(!is.character(individuals)){
            individuals <- get.individuals(x)[individuals]
        }

        ## check that all individuals are known ##
        if(!all(individuals %in% get.individuals(x))){
            temp <- paste(individuals[!individuals %in% get.individuals(x)], collapse=", ")
            warning(paste("The following individuals were not found in the data:", temp))
            individuals <- individuals[individuals %in% get.individuals(x)]
        }

        ## escape if < 2 individuals ##
        if(length(individuals)<2) return(NULL)

        ## delete obsolete edges ##
        toRemove <- which(!network.vertex.names(x@contacts) %in%  individuals) # individuals to remove
        x@contacts <- delete.vertices(x@contacts, toRemove) # remove vertices
    }

    ## DATES FROM ... ##
    if(!is.null(date.from)){
        date.from <- na.omit(date.from)

        ## subset contacts ##
        if(inherits(x@contacts, "networkDynamic")){
            x@contacts <- get.contacts(x, from=date.from, to=Inf)
        }
    } # end subset by date.from


    ## DATES TO ... ##
    if(!is.null(date.to)){
        date.to <- na.omit(date.to)

        ## subset contacts ##
        if(inherits(x@contacts, "networkDynamic")){
            x@contacts <- get.contacts(x, from=-1, to=date.to)
        }
    } # end subset by date.to

    ## do not return graphs with no edge ##
    if(nrow(suppressWarnings(as.data.frame(x))) < 1) return(NULL)

    return(x)
}) # end subset for obkContacts





####################
## obkData method ##
####################
setMethod("subset", "obkData", function(x, individuals=NULL, locus=NULL, sequences=NULL,
                                        date.from=NULL, date.to=NULL, date.format=NULL,
                                        ...){
    ## CHECK THAT REQUESTED INFO IS THERE ##
    if(is.null(x@individuals)) {
        individuals <- NULL
    }

    if(is.null(x@dna)){
        locus <- sequences <- NULL
    }


    ## SUBSET BY INDIVIDUALS ##
    if(!is.null(individuals)){
        individuals <- na.omit(individuals)

        ## handle non-character argument ##
        if(!is.character(individuals)){
            individuals <- get.individuals(x)[individuals]
        }

        ## check that all individuals are known ##
        if(!all(individuals %in% get.individuals(x))){
            temp <- paste(individuals[!individuals %in% get.individuals(x)], collapse=", ")
            warning(paste("The following individuals were not found in the data:", temp))
            individuals <- individuals[individuals %in% get.individuals(x)]
        }

        ## subset @individuals ##
        if(!is.null(x@individuals)) x@individuals <- x@individuals[individuals, ,drop=FALSE]

        ## subset @dna ##
        if(!is.null(x@dna)) x@dna <- suppressWarnings(subset(x@dna, individuals=individuals))

        ## subset @records ##
        if(!is.null(x@records)){
            for(i in 1:length(x@records)){
                x@records[[i]] <- x@records[[i]][x@records[[i]]$"individualID" %in% individuals, ,drop=FALSE]
            }
        }

        ## subset @contacts ##
        if(!is.null(x@contacts)){
            x@contacts <- suppressWarnings(subset(x@contacts, individuals=individuals))
        }
    } # end subsetting by individuals


    ## SUBSET BY LOCUS ##
    if(!is.null(locus)){
        locus <- na.omit(locus)

        ## subset @dna ##
        if(!is.null(x@dna)) x@dna <- subset(x@dna, locus=locus)

        ## keep only relevant individuals ##
        x <- suppressWarnings(subset(x, individuals=get.individuals(x@dna)))

        ## keep only relevant dates ##
        x <- subset(x, date.from=min(get.dates(x@dna),na.rm=TRUE), date.to=max(get.dates(x@dna),na.rm=TRUE))

        ## subset @contacts ##
        if(!is.null(x@contacts)){
            x@contacts <- suppressWarnings(subset(x@contacts, individuals=get.individuals(x@dna)))
        }

    } # end subsetting by locus


    ## SUBSET BY SEQUENCES ##
    if(!is.null(sequences)){
        sequences <- na.omit(sequences)

        ## subset @dna ##
        if(!is.null(x@dna)) x@dna <- subset(x@dna, sequences=sequences)

        ## keep only relevant individuals ##
        x <- suppressWarnings(subset(x, individuals=get.individuals(x@dna)))

        ## keep only relevant dates ##
        x <- subset(x, date.from=min(get.dates(x@dna),na.rm=TRUE), date.to=max(get.dates(x@dna),na.rm=TRUE))

        ## subset @contacts ##
        if(!is.null(x@contacts)){
            x@contacts <- suppressWarnings(subset(x@contacts, individuals=get.individuals(x@dna)))
        }

    } # end subsetting by sequences


    ## DATES FROM ... ##
    if(!is.null(date.from)){
        date.from <- na.omit(date.from)

        ## process date ##
        date.from <- .process.Date(date.from, format=date.format)

        ## subset @dna ##
        if(!is.null(x@dna)) x@dna <- subset(x@dna, date.from=date.from)

        ## subset records ##
        if(!is.null(x@records)){
            for(i in 1:length(x@records)){
                x@records[[i]] <- x@records[[i]][x@records[[i]]$"date" >= date.from, ,drop=FALSE]
            }
        }

        ## subset contacts ##
        if(!is.null(x@contacts) && inherits(x@contacts@contacts, "networkDynamic")){
            x@contacts <- subset(x@contacts, date.from=date.from)
        }
    } # end subset by date.from


    ## DATES TO ... ##
    if(!is.null(date.to)){
        date.to <- na.omit(date.to)

        ## process date ##
        date.to <- .process.Date(date.to, format=date.format)

        ## subset @dna ##
        if(!is.null(x@dna)) x@dna <- subset(x@dna, date.to=date.to)

        ## subset records ##
        if(!is.null(x@records)){
            for(i in 1:length(x@records)){
                x@records[[i]] <- x@records[[i]][x@records[[i]]$"date" <= date.to, ,drop=FALSE]
            }
        }

        ## subset contacts ##
        if(!is.null(x@contacts) && inherits(x@contacts@contacts, "networkDynamic")){
            x@contacts <- subset(x@contacts, date.to=date.to)
        }
    } # end subset by date.to



    ## SUBSET @TREES ##
    ## these only depend on @sample, so only one code is useful here
    if(!is.null(x@trees)){
        for(i in 1:length(x@trees)){
            toDrop <- x@trees[[i]]$tip.label[! x@trees[[i]]$tip.label %in% get.sequences(x)]
            ## check that we don't try to have a tree with one tip - creates an error
            if(length(setdiff(x@trees[[i]]$tip.label, toDrop))==1){
                x@trees <- x@trees[-i]
            } else {
                temp <- drop.tip(x@trees[[i]], tip=toDrop, trim.internal=TRUE, subtree=FALSE)
                if(is.null(temp)) {
                    x@trees <- x@trees[-i]
                } else {
                    x@trees[[i]] <- drop.tip(x@trees[[i]], tip=toDrop, trim.internal=TRUE, subtree=FALSE)
                }
            }
        }

        ## set slot to NULL if no tree left
        if(length(x@trees)==0) x@trees <- NULL
    }

    return(x)
}) # end obkData method




