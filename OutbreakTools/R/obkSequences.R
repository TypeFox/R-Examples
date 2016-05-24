
############################
####  CLASSE DEFINITION ####
############################

## CLASS DESCRIPTION:
## - instance of obkSequences store list of DNA alignments
## - sequences are stored as a (possibly named) list
## - each element of the list corresponds to a locus
## - names of the list are locus names
## - each element of the list is a DNAbin matrix
## (i.e., if there are several sequences, they are to be aligned)

setClass("obkSequences", representation(dna="listOrNULL", meta="data.frameOrNULL"),
         prototype(dna=NULL, meta=NULL))

setClassUnion("obkSequencesOrNULL", c("obkSequences", "NULL"))






######################
####  CONSTRUCTOR ####
######################

## INPUT DESCRIPTION:
## >> dna <<
## - a list of DNA sequence matrices, in DNAbin format (one matrix per locus)
## - names of the list are names of loci
## - dates and individualID are mandatory, but may be fetched from labels
## - sequence labels follow: [sequenceID][sep][individualID][sep][date]
## (note: by default, [sep] is "_")
##
## >> ... <<
## - any (possibly named) vector of meta data
## - the length and order is supposed to match sequences in '@dna'
##
setMethod("initialize", "obkSequences", function(.Object, dna=NULL, individualID=NULL,
                                                 date=NULL, ..., date.format=NULL, quiet=FALSE,
                                                 sep="_", keep.simple.labels=TRUE) {

    ## RETRIEVE PROTOTYPED OBJECT ##
    x <- .Object
    other <- NULL


    ## HANDLE DNA ##
    ## cases where an obkSequences is provided ##
    if(inherits(dna, "obkSequences")){
        individualID <- dna@meta$individualID
        date <- dna@meta$date
        if(ncol(dna@meta)>2) other <- dna@meta[,-(1:2),drop=FALSE]
        dna <- dna@dna
    }

    ## cases where no info provided ##
    if(is.null(dna)) return(x)
    if(is.matrix(dna)) dna <- list(dna)

    ## coerce items in DNA to matrices ##
    dna <- lapply(dna, as.matrix)
    NSEQ <- sum(sapply(dna, nrow))
    if(NSEQ==0){
        x@dna <- NULL
        x@meta <- NULL
        return(x)
    }

    ## convert matrices of characters into DNAbin ##
    NLOC <- length(dna)
    for(i in 1:NLOC){
        if(is.character(dna[[i]])) dna[[i]] <- as.DNAbin(dna[[i]])
    }

    ## replace with generic names if needed ##
    if(is.null(names(dna))) names(dna) <- paste("locus", 1:NLOC, sep=".")


    ## HANDLE LABELS ##
    ## extract labels ##
    labels <- unlist(lapply(dna, rownames))
    if(is.null(labels) || length(labels)!=NSEQ){
        if(!quiet) cat("\n[obkSequences constructor] missing/incomplete labels provided - using generic labels.\n")
        labels <- paste("sequence", 1:NSEQ, sep=".")

        ## assign labels ##
        temp.lab <- labels
        for(i in 1:NLOC){
            rownames(dna[[i]]) <- temp.lab[1:nrow(dna[[i]])]
            temp.lab <- temp.lab[-(1:nrow(dna[[i]]))]
        }
    }

    ## extract individualID and date if needed ##
    if(is.null(individualID) || is.null(date)){
        NFIELDS <- length(unlist(strsplit(labels[1], split=sep, fixed=TRUE)))
        if(NFIELDS>=3){
            if(!quiet) cat("\n[obkSequences constructor] extracting individualID and dates from labels.\n")

            temp <- strsplit(labels, split=sep, fixed=TRUE)
            if(!all(sapply(temp, length) == NFIELDS)) {
                warning("[obkSequences constructor] Improper labels (varying numbers of fields)")
                cat("\nCulprits are:\n")
                print(labels[sapply(temp, length) != NFIELDS])
            }
            temp <- matrix(unlist(temp), nrow=length(labels), byrow=TRUE)

            ## get information ##
            labels <- temp[,1]
            individualID <- temp[,2]
            date <- temp[,3]
            if(ncol(temp)>3) other <- as.data.frame(temp[,-(1:3),drop=FALSE])

            ## reassign labels if we need simple accession numbers ##
            if(keep.simple.labels){
                for(i in 1:NLOC){
                    rownames(dna[[i]]) <- labels[1:nrow(dna[[i]])]
                    labels <- labels[-(1:nrow(dna[[i]]))]
                }
                labels <- temp[,1]
            }
        }
    }

    ## HANDLE INDIVIDUAL ID ##
    if(is.null(individualID)){
        if(!quiet) cat("\n[obkSequences constructor] note: missing individualID\n")
        individualID <- rep(NA, NSEQ)
    }


    ## HANDLE DATE ##
    if(is.null(date)){
        if(!quiet) cat("\n[obkSequences constructor] note: missing dates\n")
        date <- rep(NA, NSEQ)
    }
    date <- .process.Date(date, format=date.format)


    ## HANDLE OTHER / ... ##
    ## retrieve information ##
    if(is.null(other)) {
        other <- list(...)
    } else {
        other <- list(other)
    }
    N.OTHER <- length(other)
    if(N.OTHER>0){
        ## use generic names if needed ##
        if(is.null(names(other))) names(other) <- paste("other", 1:N.OTHER, sep=".")

        ## convert to data.frame ##
        if(!is.data.frame(other)){
            other <- as.data.frame(other, row.names=labels)
        }

        ## check row.names, possibly reorder ##
        ## right names, possibly bad order
        if(all(labels %in% row.names(other)) || all(row.names(other) %in% labels)){
            other <- other[labels,,drop=FALSE]
        } else { ## bad names
            row.names(other) <- labels
        }
    }


    ## GET LOCUS INFO ##
    locus <- rep(names(dna), sapply(dna, nrow))

    ## RESTORE NAs (might be considered as characters "NA") ##
    individualID[individualID=="NA"] <- NA

    ## FORM FINAL OBJECT ##
    x@dna <- dna
    x@meta <- data.frame(individualID=individualID, date=date, locus=locus)
    if(N.OTHER>0) x@meta <- cbind.data.frame(x@meta, other)
    row.names(x@meta) <- labels

    return(x)
}) # end obkSequences constructor








####################
####  ACCESSORS ####
####################

################
## get.nlocus ##
################
setMethod("get.nlocus","obkSequences", function(x, ...){
    if(is.null(x@dna)) return(0)
    return(length(x@dna))
})



################
## get.locus ##
################
setMethod("get.locus","obkSequences", function(x, ...){
    if(is.null(x@dna)) return(NULL)
    return(names(x@dna))
})



###################
## get.sequences ##
###################
##  (get sequence IDs)
setMethod("get.sequences","obkSequences", function(x, ...){
    if(is.null(x)) return(NULL)
    return(unlist(lapply(x@dna, rownames)))
})



####################
## get.nsequences ##
####################
setMethod("get.nsequences","obkSequences", function(x, what=c("total","bylocus"), ...){
    what <- match.arg(what)
    nLoc <- get.nlocus(x)
    if(nLoc==0) return(0)

    temp <- sapply(x@dna, nrow)
    if(what=="bylocus") return(temp)
    return(sum(temp))
})



#####################
## get.individuals ##
#####################
setMethod("get.individuals","obkSequences", function(x, ...){
    if(is.null(x)) return(NULL)
    return(unique(x@meta$individualID))
})



######################
## get.nindividuals ##
######################
setMethod("get.nindividuals","obkSequences", function(x, ...){
    if(is.null(x)) return(0)
    return(length(get.individuals(x)))
})



###############
## get.dates ##
###############
setMethod("get.dates","obkSequences", function(x, ...){
    if(is.null(x)) return(NULL)
    return(unique(x@meta$date))
})



################
## get.ndates ##
################
setMethod("get.ndates","obkSequences", function(x, ...){
    if(is.null(x)) return(0)
    return(length(get.dates(x)))
})



#############
## get.dna ##
#############
## returns a matrix of dna sequences for a given locus
setMethod("get.dna","obkSequences", function(x, locus=NULL, id=NULL, ...){
    ## return NULL if no info ##
    nLoc <- get.nlocus(x)
    if(nLoc==0) return(NULL)

    ## RETURN SLOT CONTENT AS IS IF NOTHING ELSE ASKED ##
    if(is.null(locus) && is.null(id)) return(x@dna)

    ## INFO REQUESTED PER LOCUS ##
    if(is.null(id)){
        ## return only locus if nLoc==1 and no info on locus ##
        if(nLoc==1 && is.null(locus)) return(x@dna[[1]])

        ## otherwise use locus info ##
        if(nLoc>1 && is.null(locus)) stop("locus must be specified (data contain more than one locus)")
        return(x@dna[locus])
    }

    ## INFO REQUESTED PER SEQUENCE ID ##
    ## if logicals or integers, find corresponding names
    if(is.logical(id) | is.numeric(id) | is.integer(id)){
        id <- get.sequences(x)[id]
    }
    id <- as.character(id)
    if(!all(id[!is.na(id)] %in% get.sequences(x))) {
        temp <- paste(id[!is.na(id) & !id %in% get.sequences(x)], collapse=", ")
        warning(paste("The following sequence IDs are not in the dataset:", temp))
        id <- id[!is.na(id) & id %in% get.sequences(x)]
    }
    out <- lapply(x@dna, function(e) e[id[id %in% rownames(e)],,drop=FALSE])
    out <- out[sapply(out, nrow)>0]
    return(out)
})








######################
####  SHOW METHOD ####
######################
setMethod ("show", "obkSequences", function(object){
    nLoc <- get.nlocus(object)
    nSeq <- get.nsequences(object)
    seqword <- ifelse(nSeq>1, "sequences", "sequence")
    locword <- ifelse(nLoc>1, "loci", "locus")

    cat(paste("= @dna =\n",sep=""))
    cat(paste("[", nSeq,"DNA", seqword, "in", nLoc, locword,"]\n"))
    if(nLoc>0) print(object@dna)

    cat(paste("\n= @meta =\n",sep=""))
    cat("[ meta information on the sequences ]\n")
    if(nLoc>0) {
        if(nrow(object@meta)>10){
            print(head(object@meta,4))
            cat("\n...\n")
            print(tail(object@meta,4))
        } else {
            print(object@meta)
        }
    }
})






##################
####  TESTING ####
##################
## NOTE: THIS MUST BE COMMENTED WHEN COMPILING/INSTALLING THE PACKAGE

## library(ape)
## data(woodmouse)

## ## test constructor / show
## new("obkSequences") # empty object
## new("obkSequences", woodmouse) # no locus info
## new("obkSequences", as.matrix(woodmouse), locus=rep(c('loc1', 'loc2', 'locXX'), c(10,4,1)))


## ## test accessors
## x <- new("obkSequences", as.matrix(woodmouse), locus=rep(c('loc1', 'loc2', 'locXX'), c(10,4,1)))
## get.dna(x, locus=1)
## get.dna(x, locus="locXX")
## get.nlocus(x)
## get.nsequences(x)
