############################
####  CLASSE DEFINITION ####
############################

## CLASS DESCRIPTION:
## Instance of uniqSequences store sequences; its content includes:
## - @uniqID: data about the identical sequences, stored as a list of vectors of sequenceID
## - @uniqdna: unique dna sequences, stored as a DNAbin

## setOldClass("DNAbin")
setClass("uniqSequences", representation(uniqID="list", uniqdna="DNAbin"),prototype(uniqID=NULL, uniqdna=NULL))


######################
####  CONSTRUCTOR ####
######################

## INPUT DESCRIPTION:
## 'uniqID': a list of vectors, each vector contains the sequenceID of sequences that are identical
## 'uniqdna': a DNAbin with named unique sequences according to the list names
##
setMethod("initialize", "uniqSequences", function(.Object, uniqID=NULL, uniqdna=NULL){
    ## RETRIEVE PROTOTYPED OBJECT ##
    x <- .Object

    ## store old option ##
    o.opt <- options("stringsAsFactors")
    options("stringsAsFactors"=FALSE)
    on.exit(options(o.opt))

    ## PROCESS INFORMATION ABOUT uniqIDs ('uniqID') ##
    if (!is.null(uniqdna) && inherits(uniqdna,"DNAbin")){
      x@uniqdna<-uniqdna
    }
    if(!is.null(uniqID)){
      x@uniqID <- list()
      for(i in 1:length(uniqID))
        {
          x@uniqID[[i]] <- uniqID[[i]]
        }
      names(x@uniqID) <- names(uniqID)
      ## make sure that all the uniqIDs are in 'uniqdna'
      if(!is.null(x@uniqID) && !is.null(x@uniqdna)){
        # print(labels(x@uniqdna))
        # print(names(x@uniqID))
        unknownIDs <- unique(labels(x@uniqdna))[!unique(labels(x@uniqdna)) %in% names(x@uniqID)]
        if(length(unknownIDs)>0) {
          unknownIDs.txt <- paste(unknownIDs, collapse = ", ")
          warning(paste("the following uniqIDs in the DNAbin do not have information in the list:\n", unknownIDs.txt))
        }
        unknownlabs <- unique(labels(x@uniqID))[!unique(names(x@uniqID)) %in% labels(x@uniqdna)]
        if(length(unknownlabs)>0) {
          unknownlabs.txt <- paste(unknownlabs, collapse = ", ")
          warning(paste("the following uniqIDs in the list do not have sequences in the DNAbin:\n", unknownlabs.txt))
        }

      }
    }

    ## RETURN OBJECT ##
    return(x)
}) # end uniqSequences constructor


##################
####  TESTING ####
##################
## NOTE: THIS MUST BE COMMENTED WHEN COMPILING/INSTALLING THE PACKAGE

## new("uniqSequences",uniqID=summary.seq,uniqdna=uniqdna)
## removed an uniqID from the list
## short.summary.seq<-summary.seq
## short.summary.seq[1]<-NULL
## new("uniqSequences",uniqID=short.summary.seq,uniqdna=uniqdna)
## adding a uniqID to a list
## long.summary.seq<-summary.seq
## long.summary.seq[['uniqseqID1000']]<-paste("originalID",1:10,sep="")
## new("uniqSequences",uniqID=long.summary.seq,uniqdna=uniqdna)

## us<-new("uniqSequences",uniqID=summary.seq,uniqdna=uniqdna)
