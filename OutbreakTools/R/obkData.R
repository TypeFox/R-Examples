
############################
####  CLASSE DEFINITION ####
############################

## CLASS DESCRIPTION:
## Instance of obkData store outbreak data; its content includes:
## - @individuals: meta-information on the individuals (group, etc.), stored as a data.frame
## - @records: list of data.frames with records made on individuals; mantatory columns: individualID, date
## - @dna: dna data, stored as a obkSequences object
## - @contacts: contact information as obkContacts
## - @context: list of data.frames with records of events made on population; mantatory columns: date
setClass("obkData", representation(individuals="data.frameOrNULL", records="listOrNULL",
                                   dna="obkSequencesOrNULL", contacts="obkContactsOrNULL",
                                   context="listOrNULL", trees="multiPhyloOrNULL"),
         prototype(individuals=NULL, records=NULL, dna=NULL, contacts=NULL, trees=NULL, context=NULL))







######################
####  CONSTRUCTOR ####
######################

## INPUT DESCRIPTION:
## 'individuals': a data.frame with any information on the individuals, each row being an individual, with the following columns:
## - "individualID"
## - any other named column
##
## 'records': a list of data.frame where each row is an observation made on an individual, and the following mandatory columns:
## - "individualID"
## - "date"
##
## 'dna': a DNAbin list with named sequences
##
## 'records': list of clinical datasets, each stored as a data.frame
##
## 'contacts': a matrix of characters indicating edges using two columns; if contacts are directed,
## the first column is 'from', the second is 'to'; values should match individual IDs (as returned
## by get.individuals(x)); if numeric values are provided, these are converted as integers and
## assumed to correspond to individuals returned by get.individuals(x).
##
## 'contacts.start': a vector of dates indicating the beginning of each contact
##
## 'contacts.end': a vector of dates indicating the end of each contact
##
## 'contacts.duration': another way to specify contacts.end, as duration of contact
##
## 'context': a list of data.frame where each row is an observation made on a population (which includes interventions,
## environmental observations etc.), and the following mandatory columns:
## - "date"
##
setMethod("initialize", "obkData", function(.Object, individuals=NULL, records=NULL, dna=NULL,
                                            trees=NULL, contacts=NULL, context=NULL,
                                            contacts.start=NULL, contacts.end=NULL,
                                            contacts.duration=NULL, contacts.directed=FALSE,
                                            date.format=NULL, dna.individualID=NULL, dna.date=NULL,
                                            dna.date.format=date.format, dna.sep="_", quiet=FALSE,
                                            check=TRUE, ...){

    ## RETRIEVE PROTOTYPED OBJECT ##
    x <- .Object


    ## SORT OUT STRINGSASFACTORS  OPTION ##
    o.opt <- options("stringsAsFactors")
    options("stringsAsFactors"=FALSE)
    on.exit(options(o.opt))


    ## HANDLE OBKDATA INPUT ##
    if(inherits(individuals, "obkData")){
        individuals <- x@individuals
        records <- x@records
        dna <- x@dna
        trees <- x@trees
        contacts <- x@contacts
        context <- x@context
    }


    ## HANDLE INDIVIDUALS ##
    ## force NULL if empty data.frame ##
    if(!is.null(individuals) && (nrow(individuals)==0 || ncol(individuals)==0)) individuals <- NULL

    ## process information ##
    if(!is.null(individuals)) {
        ## force type to data.frame
        individuals <- as.data.frame(individuals)

        ## check mandatory fields
        if(!"individualID" %in% names(individuals) && is.null(row.names(individuals)))
            stop("no field 'individualID' in the individuals data.frame ('individuals')")

        ## extract/process labels
        if("individualID" %in% names(individuals)){
            lab <- as.character(individuals[,"individualID"])
        } else {
            lab <- as.character(row.names(individuals))
        }

        ## store info in output
        x@individuals <- individuals[, names(individuals)!="individualID", drop=FALSE]
        row.names(x@individuals) <- lab
    } else {
        x@individuals <- NULL
    }


    ## HANDLE RECORDS ##
    ## force NULL if empty list ##
    if(!is.null(records) && length(records)==0) records <- NULL

    ## force NULL if empty data.frame ##
    if(!is.null(records) && is.data.frame(records) && (nrow(records)==0 || ncol(records)==0)) records <- NULL

    ## process information ##
    if(!is.null(records)) {
        ## force type to list
        if(is.data.frame(records)) records <- list(records)

        ## remove NULL elements
        records <- records[!sapply(records, is.null)]

        ## check mandatory fields
        NREC <- length(records)
        for(i in 1:NREC){
            if(!"individualID" %in% names(records[[i]])) stop(paste("no field 'individualID' in the records data.frame", names(records)[i], ")"))
            if(!"date" %in% names(records[[i]])) stop(paste("no field 'date' in the records data.frame", names(records)[i], ")"))
        }

        ## store info in output
        ## (reorder the columns within each data frame / convert types)
        all.records.ID <- NULL
        for(i in 1:NREC){
            nameOrder <- c(c("individualID","date"), setdiff(names(records[[i]]), c("individualID","date")))
            x@records[[i]] <- records[[i]][, nameOrder]
            x@records[[i]][,"individualID"] <- as.character(x@records[[i]][,"individualID"])
            if(is.factor(x@records[[i]][,"date"])) x@records[[i]][,"date"] <- as.character(x@records[[i]][,"date"])
            x@records[[i]][,"date"] <- .process.Date(x@records[[i]][,"date"], format=date.format)
            all.records.ID <- c(all.records.ID, x@records[[i]][, "individualID"])
        }

        names(x@records) <- names(records)
    } else { # no information
        x@records <- NULL
    }


    ## HANDLE DNA ##
    if(!is.null(dna)){
        ## if dna is already an obkSequences
        if(inherits(dna, "obkSequences")){
            x@dna <- dna
        } else {
            ## pass on inputs to obkSequences constructor
            x@dna <- new("obkSequences",
                         dna=dna, individualID=dna.individualID, date=dna.date,
                         ..., date.format=dna.date.format, quiet=quiet, sep=dna.sep)
        }
    }


    ## HANDLE PHYLOGENIES ('trees') ##
    if(!is.null(trees)){
        ## check class
        if(!inherits(trees, "multiPhylo")) stop("trees must be a multiPhylo object")

        ## check label consistency (to be added)
        x@trees <- trees
    }


    ## HANDLE CONTACTS ##
    if(!is.null(contacts)){
        if(inherits(contacts, "obkContacts")){
            x@contacts <- contacts
        } else {
            ## process vertices provided as numbers
            if(is.numeric(contacts)){
                if(!is.null(x@individuals)){
                    ## replace with labels if available
                    contacts <- matrix(row.names(x@individuals)[contacts], ncol=2)
                } else {
                    ## convert as characters otherwise
                    contacts <- matrix(as.character(contacts), ncol=2)
                }
            }

            ## pass arguments to the obkContacts constructor
            x@contacts <- new("obkContacts", from=contacts[,1,drop=TRUE], to=contacts[,2,drop=TRUE],
                              directed=contacts.directed, start=contacts.start, end=contacts.end,
                              duration=contacts.duration)
        }
    }


    ## HANDLE CONTEXT ##
    ## force NULL if empty list ##
    if(!is.null(context) && length(context)==0) context <- NULL

    ## force NULL if empty data.frame ##
    if(!is.null(context) && is.data.frame(context) && (nrow(context)==0 || ncol(context)==0)) context <- NULL

    ## process information ##
    if(!is.null(context)) {
      ## force type to list
      if(is.data.frame(context)) context <- list(context)

      ## remove NULL elements
      context <- context[!sapply(context, is.null)]

      ## check mandatory fields
      NCONT <- length(context)
      for(i in 1:NCONT){
        if(!"date" %in% names(context[[i]])) stop(paste("no field 'date' in the context data.frame", names(context)[i], ")"))
      }

      ## store info in output
      ## (reorder the columns within each data frame / convert types)
      for(i in 1:NCONT){
        nameOrder <- c(c("date"), setdiff(names(context[[i]]), c("date")))
        x@context[[i]] <- context[[i]][, nameOrder]
        if(is.factor(x@context[[i]][,"date"])) x@context[[i]][,"date"] <- as.character(x@context[[i]][,"date"])
        x@context[[i]][,"date"] <- .process.Date(x@context[[i]][,"date"], format=date.format)
      }

      names(x@context) <- names(context)
    } else { # no information
      x@context <- NULL
    }


    ## QUALITY/CONSISTENCY CHECKS ##
    if(check){
        ## look for undocumented individuals in @records ##
        if(!is.null(x@individuals) && !is.null(x@records)){
            temp <- all.records.ID[!is.na(all.records.ID)]
            unknownIDs <- unique(temp)[!unique(temp) %in% row.names(x@individuals)]
            if(length(unknownIDs)>0) {
                unknownIDs.txt <- paste(unknownIDs, collapse = ", ")
                warning(paste("records refer to undocumented individuals:\n", unknownIDs.txt))
            }
        }

        ## look for undocumented individuals in @dna ##
        if(!is.null(x@individuals) && !is.null(x@dna)){
            temp <- x@dna@meta$individualID[!is.na(x@dna@meta$individualID)]
            unknownIDs <- temp[!temp %in% row.names(x@individuals)]
            if(length(unknownIDs)>0){
                unknownIDs.txt <- paste(unknownIDs, collapse = ", ")
                warning(paste("dna sequences refer to undocumented individuals:\n", unknownIDs.txt))
            }
        }
    }


    ## RETURN OBJECT ##
    return(x)
}) # end obkData constructor











##################
####  TESTING ####
##################
## NOTE: THIS MUST BE COMMENTED WHEN COMPILING/INSTALLING THE PACKAGE

## ## EMPTY OBJECT ##
## new("obkData")


## ## INDIVIDUAL INFO ONLY ##
## new("obkData", individuals=data.frame("individualID"=letters))
## new("obkData", individuals=data.frame("individualID"=letters, age=1:26, 1:26))

## ## INDIVIDUAL INFO ONLY ##
## new("obkData", individuals=data.frame("individualID"=letters))
## new("obkData", individuals=data.frame("individualID"=letters, age=1:26, 1:26))


## ## DNA INFO, NOTHING ELSE ##
## library(ape)
## data(woodmouse)
## dat.dna <- as.list(woodmouse)

## new("obkData", dna=dat.dna) # should be empty

## ## TEST CONTEXT

# library(OutbreakTools)
# data(ToyOutbreakRaw)
# names(ToyOutbreakRaw)

# vaccination <- data.frame(date=c("02/01/2000","04/01/2000"),Targeted_population=c("children","elderly"))
# school_closure <-data.frame(date=c("05/01/2000"),ending_date=c("26/01/2000"), location=c("London"))
# temperature <-data.frame(date=c("02/01/2000","03/01/2000","04/01/2000","05/01/2000","06/01/2000","07/01/2000"),temp=c(4,2,-1,2,0,1))
#
# context <-list(vaccination,school_closure,temperature)
# names(context) <-c("vaccination","school_closure","temperature")
#
# x <- new("obkData", individuals = ToyOutbreakRaw$individuals, context=context,date.format = "%d/%m/%Y")

## and test accessors
# get.context(x)
# get.ncontext(x)
# get.data(x,"context")
# get.data(x,"date",where="context")
# get.data(x,"temperature",where="context")
# get.data(x,"vaccination",where="context")
# get.data(x,"vaccination")
# get.data(x,"date")
# get.data(x,"temperature")
