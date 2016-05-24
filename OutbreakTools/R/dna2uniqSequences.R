## Function to convert a class 'DNAbin' to the class 'uniqSequences'
##
## param dna an object of the class "DNAbin"
## author Joseph Hughes
dna2uniqSequences<-function(dna){
    if(is.null(dna)) return(NULL)
    if(!inherits(dna,"DNAbin")) {
        warning("dna should be a DNAbin object - returning NULL")
        return(NULL)
    }
    seqmatrix <- as.character(dna)
    seqstring<-apply(format(seqmatrix), 1, paste, collapse="")
    uniqList<-tapply(names(seqstring),list(seqstring),I)
    splitseqstr<-strsplit(names(uniqList),"")
    uniqmat<-do.call(rbind, splitseqstr)
    uniqnames<-paste("uniqseqID",seq(1:length(uniqList)),sep="")
    rownames(uniqmat)<-uniqnames
    names(uniqList)<-uniqnames
    uniqdna<-as.DNAbin(uniqmat)
    uniqSequences<-new("uniqSequences",uniqID=uniqList,uniqdna=uniqdna)


  return(uniqSequences)
}

##################
####  TESTING ####
##################
## NOTE: THIS MUST BE COMMENTED WHEN COMPILING/INSTALLING THE PACKAGE

## Might want to write a get.dnasubset accessor for this
## get uniq sequences
## Extracting sequenceIDs from the samples table
## seqids<-na.omit(samples[samples$sampleID=="904","sequenceID"])
## get the index of the sequenceIDs in the DNAbin
## seqindex<-which(labels(dna) %in% seqids)
## create a subset DNAbin
## subsetdna<-dna[seqindex, ]
## get a particular sequence id with summary.seq$uniqseqID4[5]
## get number of sequences length(summary.seq$uniqseqID4)
## uniqsubset<-dna2uniqSequences(subsetdna)
