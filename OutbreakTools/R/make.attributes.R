##############################
## make.sequence.attributes ##
##############################
make.sequence.attributes <- function(x){
    ## check object type##
    if(!inherits(x,"obkData")) stop("x is not an obkData object")

    ## start with what is in @dna@meta ##
    out <- x@dna@meta

    ## important to keep track of the labels as merge shuffles rows around ##
    out$sequenceID <- rownames(out)

    ## merge with info on individuals (@individuals) ##
    if(!is.null(x@individuals)){
        df.individuals <- get.data(x, "individuals")
        df.individuals$individualID <- rownames(df.individuals)
        out <- merge(out, df.individuals, by = "individualID", all.x = TRUE, sort=FALSE)
    }

    ## merge with info on records (@records) ##
    if(!is.null(x@records)){
        ## merge each data.frame (by individualID)
        for(i in 1:length(x@records)){
            temp <- x@records[[i]]
            ##names(temp) <- gsub("date", paste(names(x@records)[i],"date", sep="."), names(temp))
            out <- merge(out, temp, all.x = TRUE, sort=FALSE)
        }
    }


    ## restore original row.names and ordering ##
    rownames(out) <- out$sequenceID
    out$sequenceID <- NULL
    out <- out[rownames(x@dna@meta),,drop=FALSE]

    return(out)
} # end make.sequence.attributes






#########################
## make.tip.attributes ##
#########################
make.tip.attributes <- function(x, which.tree = 1){
    ## check object type##
    if(!inherits(x,"obkData")) stop("x is not an obkData object")

    ## start with what is in @dna@meta ##
    phylo <- get.trees(x)[[which.tree]]
    out <- x@dna@meta[phylo$tip.label, , drop=FALSE]

    ## important to keep track of the labels as merge shuffles rows around ##
    out$tip.label <- rownames(out)

    ## merge with info on individuals (@individuals) ##
    if(!is.null(x@individuals)){
        df.individuals <- get.data(x, "individuals")
        df.individuals$individualID <- rownames(df.individuals)
        out <- merge(out, df.individuals, by = "individualID", all.x = TRUE, all.y=FALSE, sort=FALSE)
    }

    ## merge with info on records (@records) ##
    if(!is.null(x@records)){
        ## merge each data.frame (by individualID)
        for(i in 1:length(x@records)){
            temp <- x@records[[i]]
            ##names(temp) <- gsub("date", paste(names(x@records)[i],"date", sep="."), names(temp))
            out <- merge(out, temp, all.x = TRUE, all.y=FALSE, sort=FALSE)
        }
    }


    ## restore original row.names and ordering ##
    rownames(out) <- out$tip.label
    out$tip.label <- NULL
    out <- out[phylo$tip.label,,drop=FALSE]

    return(out)
} # end make.tip.attributes








################################
## make.individual.attributes ##
################################
make.individual.attributes <- function(x){
    ## check object type##
    if(!inherits(x,"obkData")) stop("x is not an obkData object")

    ## start with what is in @dna@meta ##
    out <- x@individuals

    ## important to keep track of the labels as merge shuffles rows around ##
    out$individual.label <- out$individualID <- rownames(out)

    ## merge with info on records (@records) ##
    if(!is.null(x@records)){
        ## merge each data.frame (by individualID)
        for(i in 1:length(x@records)){
            temp <- x@records[[i]]

            ## need to reshape data.frame into wide format when multiple individuals ##
            if(any(table(temp$individualID)>1)){
                newtab <- split(temp, temp$individualID)
                for(j in 1:length(newtab)) newtab[[j]]$time <- order(newtab[[j]]$date)
                newtab <- Reduce(rbind.data.frame, newtab)
                temp <- reshape(newtab, idvar="individualID", direction="wide")
            }

            ## alter the names ##
            toChange <- names(temp)!="individualID"
            names(temp)[toChange] <- paste(names(x@records)[i], names(temp)[toChange], sep=".")

            ## merge ##
            out <- merge(out, temp, all.x = TRUE, all.y=FALSE, by="individualID", sort=FALSE)
        }
    }

     ## merge with info on dna (@dna@meta) ##
    if(!is.null(x@dna)){

        ## keep only unique indiv/date combinations ##
        splitDates <- tapply(as.character(x@dna@meta$date), x@dna@meta$individualID, unique)
        temp <- data.frame(date=.process.Date(unlist(splitDates)))
        temp$individualID <- rep(names(splitDates), sapply(splitDates, length))

        ## need to reshape data.frame into wide format when multiple individuals ##
        if(any(table(temp$individualID)>1)){
            newtab <- split(temp, temp$individualID)
            for(j in 1:length(newtab)) newtab[[j]]$time <- order(newtab[[j]]$date)
            newtab <- Reduce(rbind.data.frame, newtab)
            temp <- reshape(newtab, idvar="individualID", direction="wide")
        }

        ## alter the names ##
        toChange <- names(temp)!="individualID"
        names(temp)[toChange] <- paste("dna", names(temp)[toChange], sep=".")

        ## merge ##
        out <- merge(out, temp, all.x = TRUE, all.y=FALSE, by="individualID", sort=FALSE)
    }


    ## restore original row.names and ordering ##
    rownames(out) <- out$individual.label
    out$individual.label <- NULL
    out <- out[rownames(x@individuals),,drop=FALSE]

    return(out)
} # end make.individual.attributes
