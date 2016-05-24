
########################
####  BASIC METHODS ####
########################


##########
## show ##
##########
setMethod("show", "obkData", function(object){
    N <- length(slotNames(object))
    cat("\n=== obkData object ===")
    empty <- rep(TRUE, N)
    for(i in 1:N){
        if(!is.null(slot(object, slotNames(object)[i]))){
            cat(paste("\n== @", slotNames(object)[i], " == \n",sep=""))
            print(slot(object, slotNames(object)[i]))
            empty[i] <- FALSE
        }
    }

    if(any(empty)){
        txt <- paste("@", slotNames(object)[empty], collapse=", ", sep="")
        cat("\n== Empty slots == \n", txt)
    }

    cat("\n")
})




#############
## summary ##
#############
setMethod("summary", "obkData", function(object, ...){
    ## FUNCTION TO DISPLAY SUMMARY OF 1 DATA.FRAME ##
    ## displays the numer of entries, individuals, and time window for a given data.frame
    ## plus summary of other variables (not individualID and date)
    f1 <- function(x, indent="  ", ...){
        if(is.null(x$individualID)){
            cat(indent, nrow(x), " entries\n", sep="")
        } else {
          if(!all(is.na(x$date))){
            cat(indent, nrow(x), " entries,  ", length(unique(x$individualID)), " individuals, from ",
                as.character(min(x$date, na.rm=TRUE)), " to ", as.character(max(x$date,na.rm=TRUE)), "\n", sep="")
          }else{
            cat(indent, nrow(x), " entries,  ", length(unique(x$individualID)), " individuals, from ",
                as.character(min(x$date)), " to ", as.character(max(x$date)), "\n", sep="")
          }
        }
        if(ncol(x)>2){
            temp <- x
            temp$individualID <- NULL
            temp$date <- NULL
            cat(indent, "recorded fields are:\n", sep="")
            for(i in 1:ncol(temp)){
                cat(indent, "<", names(temp)[i], "> ", sep="")
                .inlineSummary(temp[[i]], ...) #  '...' is passed to 'format'
            }
        }
    }


    cat(paste("Dataset of ",get.nindividuals(object,"all")," individuals with...\n",sep=""))

    ## handle @individuals ##
    if(!is.null(object@individuals)){
        cat("== @individuals ==\n")
        cat("individuals information\n")
        f1(object@individuals)
        cat("\n")
    }

    ## handle @records ##
    if(!is.null(object@records)){
        temp <- paste(get.records(object), collapse=", ")
        cat("== @records ==\n")
        cat("records on: ", temp,"\n")
        for(i in 1:length(object@records)){
            cat("$", names(object@records)[i], "\n", sep="")
            f1(object@records[[i]])
        }
        cat("\n")
    }

    ## handle @context ##
    if(!is.null(object@context)){
        temp <- paste(get.context(object), collapse=", ")
        cat("== @context ==\n")
        cat("contextual information including", temp,"\n")
        for(i in 1:length(object@context)){
            cat("$", names(object@context)[i], "\n", sep="")
            f1(object@context[[i]])
        }
        cat("\n")
    }

    ## handle @dna ##
    if(!is.null(object@dna)){
        cat("== @dna ==\n")
        if(!all(is.na(object@dna@meta$date))){
        cat(get.nsequences(object)," sequences across ", get.nlocus(object), " loci, ",
            get.nindividuals(object@dna), " individuals, from ", as.character(min(object@dna@meta$date, na.rm=TRUE)),
            " to ", as.character(max(object@dna@meta$date, na.rm=TRUE)), "\n", sep="")
        }else{
          cat(get.nsequences(object)," sequences across ", get.nlocus(object), " loci, ",
              get.nindividuals(object@dna), " individuals, from ", as.character(min(object@dna@meta$date)),
              " to ", as.character(max(object@dna@meta$date)), "\n", sep="")
        }
        cat("length of concatenated alignment: ", sum(sapply(object@dna@dna,ncol)), " nucleotides\n", sep="")
        cat("Attached meta data:\n")
        f1(object@dna@meta)
        cat("\n")
    }

    ## handle @contacts ##
    if(!is.null(object@contacts)){
        cat("== @contacts ==\n")
        cat(get.ncontacts(object@contacts)," contacts between ", get.nindividuals(object@contacts), " individuals\n", sep="")
        cat("\n")
    }

    ## handle @trees ##
    if(!is.null(object@contacts)){
        cat("== @trees ==\n")
        cat(length(object@trees)," phylogenetic trees with ", length(object@trees[[1]]$tip.label), " tips\n", sep="")
        cat("\n")
    }

    return(invisible())
}) # end summary


## test:
## library(OutbreakTools)
## data(HorseFlu)
## summary(HorseFlu)
## summary(new("obkData"))






##########
## head ##
##########
setMethod("head", "obkData", function(x, n=4L, ...){
    Nslots <- length(slotNames(x))
    cat("\n=== obkData x ===")
    empty <- rep(TRUE, Nslots)
    for(i in 1:Nslots){
        if(!is.null(slot(x, slotNames(x)[i]))){ # if slot is not NULL
            cat(paste("\n== @", slotNames(x)[i], "== \n",sep=""))
            if(inherits(slot(x, slotNames(x)[i]), c("obkSequences","obkContacts","multiPhylo"))){ # special classes
                print(slot(x, slotNames(x)[i]))
            } else if(is.list(slot(x, slotNames(x)[i])) && !is.data.frame(slot(x, slotNames(x)[i]))){ # use custom 'head' for lists
                lapply(slot(x, slotNames(x)[i]), function(e) print(head(e, n=n, ...)))
            } else {
                print(head(slot(x, slotNames(x)[i]), n=n, ...))
            }
            empty[i] <- FALSE
        }
    }

    if(any(empty)){
        txt <- paste("@", slotNames(x)[empty], collapse=", ", sep="")
        cat("\n== Empty slots == \n", txt)
    }

    cat("\n")
})





##########
## tail ##
##########
setMethod("tail", "obkData", function(x, n=4L, ...){
    Nslots <- length(slotNames(x))
    cat("\n=== obkData x ===")
    empty <- rep(TRUE, Nslots)
    for(i in 1:Nslots){
        if(!is.null(slot(x, slotNames(x)[i]))){
            cat(paste("\n== @", slotNames(x)[i], "== \n",sep=""))
            if(inherits(slot(x, slotNames(x)[i]), c("obkSequences","obkContacts","multiPhylo"))){ # special classes
                print(slot(x, slotNames(x)[i]))
            } else if(is.list(slot(x, slotNames(x)[i])) && !is.data.frame(slot(x, slotNames(x)[i]))){ # use custom 'tail' for lists
                lapply(slot(x, slotNames(x)[i]), function(e) print(tail(e, n=n, ...)))
            } else {
                print(tail(slot(x, slotNames(x)[i]), n=n, ...))
            }
            empty[i] <- FALSE
        }
    }

    if(any(empty)){
        txt <- paste("@", slotNames(x)[empty], collapse=", ", sep="")
        cat("\n== Empty slots == \n", txt)
    }

    cat("\n")
})










