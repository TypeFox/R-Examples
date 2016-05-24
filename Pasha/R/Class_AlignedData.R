### Class definition

setClass("AlignedData", representation=representation(
                readID="integer",
                seqnames="factor",
                strand="factor",
                position="integer",
                mapq="integer",
                isize="integer",
                qwidth="integer",
                flag="integer",
                mrnm="factor",
                notAligned="integer",
                pairedEnds="logical",
                weight="numeric"),
#                            prototype=prototype(strand=factor(levels=.STRAND_LEVELS),
#                                                mapq=NumericQuality()),

        validity=function(object)
        {
            if(length(readID(object))>0)
            {
                # Forbid the specification of weight for pairedEnds, could be confusing in further processing, this situation has to be dealt with before (user handles the pairs and weights relation by defining subobjects)
                if((length(weight(object))>0) && pairedEnds(object)) return("AlignedData objects assume that the user deal with the relation between pairs and read weights, it is thus not possible to define weights for a dataset declared as paired-ends")
                if((length(isize(object))==0) && pairedEnds(object)) return("AlignedData objects assume that for paired-ends experiments, the insert size is specified while creating object")
                
                # Function used in sapply to get the slots length
                getSlotsLength <- function(slotName, object){return(length(slot(object, slotName)))}
                
                ### Check that all members (except atomic ones and data.frame) have the same length
                
                # These are atomic slots, designed to get stats on the whole experiment/lane
                atomicSlots <- c("notAligned", "pairedEnds")
                if(!all(sapply(atomicSlots,getSlotsLength,object)==1)) return(paste("Some elements are defined as atomic and are mandatory, please be sure to assign a unique value : ",paste(atomicSlots, collapse=", "),sep=""))
                
                
                # These are optional slots. Briefly, they will be defined for BAM and/or Paired-Ends datasets (either size=0 or same size as others) and/or reads that aligned in several position
                optionalSlots<- c("mapq", "isize", "flag", "mrnm", "weight", "qwidth")
                optionalSlotsLength <- sapply(optionalSlots, getSlotsLength,object)
                optionalSlots <- optionalSlots[optionalSlotsLength==0] # Keep as optional the ones which have length 0, the other ones needs to be checked
                
                #print(optionalSlots[optionalSlotsLength!=0])
                
                # Select 'non atomic' and 'not of length 0 and optional' members to check for their size consistency
                otherSlots <- slotNames(object) 
                otherSlots <- otherSlots[!(otherSlots %in% atomicSlots)]
                otherSlots <- otherSlots[!(otherSlots %in% optionalSlots)]
                
                slotLengths <- sapply(otherSlots, getSlotsLength,object)
                if(any(slotLengths!=length(object@readID))) return(paste("All defined elements of an AlignedData object must have the same length except atomic ones (",paste(atomicSlots, collapse=", "),")\n", sep=""))
                
                # Restore the eventually levels definition for strand to bioconductor standard...
                object@strand <- factor(object@strand, levels=c("+","-","*"))
                
                # Detect the unsupported strand levels (can be dropped using 'dropUndefinedStrand')
                undefinedStrand <- is.na(object@strand) | object@strand=="*"
                if(any(undefinedStrand))
                {
                    warning(paste("Some reads (", sum(undefinedStrand),") contain unsupported strand value (supported values : '+', '-')...", sep=""))
                }
                
                # Check that paired-end datasets have required slots and have reads marked as paired (by flag) and have isize information 
                if(pairedEnds(object))
                {
                    if(!any(as.logical(bitAnd(flag(object), 1)))) return("The object of class 'AlignedData' is declared to contain paired-end reads but not a single read has the corresponding flag bit set...")
                    if(all(is.na(isize(object)))) return("The object of class 'AlignedData' is declared to contain paired-end reads but no 'insert size' information is available (all values are 'NA')...")
                }
                
                return(TRUE)
            }
        })

        
        
### Initiator

setMethod(f="initialize", 
        signature="AlignedData",
        definition=function(.Object, readID=integer(0), seqnames=factor(), strand=factor(), position=integer(0), mapq=integer(0), isize=integer(0), qwidth=integer(0), flag=integer(0), mrnm=factor(), notAligned=integer(0), pairedEnds=logical(0), weight=double(0), drop=FALSE)
        {
            
            # afecting the values
            .Object@readID <- readID
            .Object@seqnames <- if(drop) factor(seqnames) else seqnames # Reduce the factors if necessary and requested
            .Object@strand <- strand # Reduce the factors if necessary and requested (removed since it has to be '-' '+' '*')
            .Object@position <- position
            .Object@mapq <- mapq
            .Object@isize <- isize
            .Object@qwidth <- qwidth
            .Object@flag <- flag
            .Object@mrnm <- if(drop) factor(mrnm) else mrnm # Reduce the factors if necessary and requested
            .Object@notAligned <- notAligned
            .Object@pairedEnds <- pairedEnds
            .Object@weight <- weight
            
            validObject(.Object)
            
            return(.Object)
        })

### Getters/Setters/Subsetters
# Generic function that automatically defines getters
#        setGenericSetters=function(className) # This is probably not really compliant with OOP 'good behaviour'...
#        {
#            for(currentSlot in slotNames(className))
#            {
#                # Set the generic method
#                setGeneric(name=eval(currentSlot), def=function(object, ...) {standardGeneric(currentSlot)}) # WTF !!! The 'standardGeneric' function takes the name of the argument instead of its content !!!
#                
#                # Set the method definition
#                setMethod(f=currentSlot,
#                            signature=className, 
#                            definition=function(object, ...) {return(slot(object, currentSlot))})
#            }
#        }
#        
#        setGenericSetters("AlignedData")



### Conversion to/from GAlignments objects

setAs(from="AlignedData",
        to="GAlignments",
        def=function(from) 
        {
            return(GAlignments(seqnames=seqnames(from), pos=position(from), cigar=paste(qwidth(from),"M", sep=""), strand=Rle(factor(strand(from), levels=c("+","-","*"))), names=as.character(readID(from))))
        })

setAs(from="GAlignments",
        to="AlignedData",
        def=function(from) 
        {
            # Check if there is non-standard CIGAR (insertions/deletions) and drop them with a warning
            nonStandard <- grepl("[a-ln-z]", cigar(from), ignore.case=TRUE, perl=TRUE)
            if(sum(nonStandard)>=length(from)) stop("All items from GAlignments object were reads with non-exact alignments while AlignedData objects can only contain exact matches...")
            if(any(nonStandard))
            {
                warning(paste("Some items (", sum(nonStandard),") from GAlignments object were containing non-exact alignments and were dropped to fit AlignedData standard...", sep=""))
                from <- from[!nonStandard]
            }
            # Create the AlignedData object
            return(new("AlignedData", readID=1:length(from), seqnames=as.character(seqnames(from)), strand=as.character(strand(from)), position=start(from), mapq=integer(0), qwidth=width(from), notAligned=as.integer(0), pairedEnds=FALSE))
        })

 

### Show method

# Called by show or directly by the pipeline that needs no intro (class name)
.summarizeAlignedData=function(object, noWarning=FALSE)
{
    cat("\n", length(object), " ", ifelse(length(object) == 1L, "alignment ", "alignments "), sep="")
    if(pairedEnds(object))
    {
        cat("declared as paired-ends")
        
        # Insert size distribution
        cat("\nInserts size : ")
        insertsSize <- isize(object)
        naSize <- is.na(insertsSize)                
        if(any(naSize) && (!noWarning)) warning("Some aligned reads in AlignedData object did not have any insert size declared, orphan reads and pairs consistency should be checked properly (see 'getOrphansIndexes' and 'checkPairsOK' methods)...")
        insertsSize <- insertsSize[(insertsSize>0) & (!naSize)]
        if(length(insertsSize)>0)
        {
            if(length(table(insertsSize))>1)
            {
                cat("variable...")
                insertsSizeDistribution <- summary(cut(insertsSize, breaks=5))
                resNULL <- mapply(function(intervalName, intervalCount){cat("\n", intervalName,"->", intervalCount)}, names(insertsSizeDistribution), insertsSizeDistribution)
            }
            else
            {
                cat(insertsSize[1])
            }
        }
        else
        {
            cat("Unavailable...")
        }
    }
    
    # Summarize reads size distribution
    cat("\nReads size : ")
    if(length(qwidth(object))>0)
    {
        if(length(table(qwidth(object)))>1)
        {
            cat("variable...")
            readsSizeDistribution <- summary(cut(qwidth(object), breaks=5))
            resNULL <- mapply(function(intervalName, intervalCount){cat("\n", intervalName,"->", intervalCount)}, names(readsSizeDistribution), readsSizeDistribution)
        }
        else
        {
            cat(qwidth(object)[1])
        }
    }
    else
    {
        cat("Unavailable...")
    }
    
    # Summarize alignment quality scores distribution
    cat("\nAlignment quality scores : ")
    if(length(table(mapq(object)))>1)
    {
        mapqDistribution <- summary(cut(mapq(object), breaks=5))
        resNULL <- mapply(function(intervalName, intervalCount){cat("\n", intervalName,"->", intervalCount)}, names(mapqDistribution), mapqDistribution)
    }
    else
    {
        cat("Unavailable...")
    }
    
    # Not aligned
    cat("\n", notAligned(object), " reads were reported as not aligned to reference", sep="")
    
    # Chromosomes concerned
    cat("\nChromosomes (seqnames) covered by at least a read : ", paste(unname(unique(seqnames(object))), collapse=" - "), "\n", sep="")
}

setMethod(f="show", signature="AlignedData", definition=function(object)
        {
            # Intro with class name
            cat(class(object), " object", sep="")
            # Content information
            .summarizeAlignedData(object, noWarning=FALSE)
        })



### Getters

#setGeneric(name="length") # the generic already exists in base...
setMethod(f="length",
        signature="AlignedData", 
        definition=function(x) {return(length(x@readID))})


setGeneric(name="readID", def=function(object, ...) standardGeneric("readID"))
setMethod(f="readID",
        signature="AlignedData", 
        definition=function(object, ...) {return(slot(object, "readID"))})

#setGeneric(name="seqnames") # defined in GenomeInfoDb
setMethod(f="seqnames",
        signature="AlignedData",
        definition=function(x) {return(slot(x, "seqnames"))})

#setGeneric(name="strand") # defined in BiocGenerics
setMethod(f="strand",
        signature="AlignedData", 
        definition=function(x, ...) {return(slot(x, "strand"))})

#setGeneric(name="position") # defined in ShortRead
setMethod(f="position",
        signature="AlignedData", 
        definition=function(object, ...) {return(slot(object, "position"))})

setGeneric(name="mapq", def=function(object, ...) standardGeneric("mapq"))
setMethod(f="mapq",
        signature="AlignedData", 
        definition=function(object, ...) {return(slot(object, "mapq"))})

setGeneric(name="isize", def=function(object, ...) standardGeneric("isize"))
setMethod(f="isize",
        signature="AlignedData", 
        definition=function(object, ...) {return(slot(object, "isize"))})

#setGeneric(name="qwidth") # defined in GenomicAlignments
setMethod(f="qwidth",
        signature="AlignedData", 
        definition=function(x) {return(slot(x, "qwidth"))})

#setGeneric(name="flag") # defined in ShortRead
setMethod(f="flag",
        signature="AlignedData", 
        definition=function(object, ...) {return(slot(object, "flag"))})

setGeneric(name="mrnm", def=function(object, ...) standardGeneric("mrnm"))
setMethod(f="mrnm",
        signature="AlignedData", 
        definition=function(object, ...) {return(slot(object, "mrnm"))})

setGeneric(name="notAligned", def=function(object, ...) standardGeneric("notAligned"))
setMethod(f="notAligned",
        signature="AlignedData", 
        definition=function(object, ...) {return(slot(object, "notAligned"))})

setGeneric(name="pairedEnds", def=function(object, ...) standardGeneric("pairedEnds"))
setMethod(f="pairedEnds",
        signature="AlignedData", 
        definition=function(object, ...) {return(slot(object, "pairedEnds"))})

setGeneric(name="weight", def=function(object, ...) standardGeneric("weight"))
setMethod(f="weight",
        signature="AlignedData", 
        definition=function(object, ...) {return(slot(object, "weight"))})



### Setters (The only real interest of all this is to call the validObject function after altering the data to check the object consistency before replacing the original one with the altered copy)

setGeneric(name="readID<-", def=function(object, value) standardGeneric("readID<-"))
setReplaceMethod(f="readID",
        signature="AlignedData", 
        definition=function(object, value) 
        {
            object@readID<-value
            validObject(object)
            return(object)
        })

#setGeneric(name="seqnames<-") # defined in GenomeInfoDb
setReplaceMethod(f="seqnames",
        signature="AlignedData", 
        definition=function(x, value) 
        {
            x@seqnames<-value
            validObject(x)
            return(x)
        })

#setGeneric(name="strand<-") # defined in BiocGenerics
setReplaceMethod(f="strand",
        signature="AlignedData", 
        definition=function(x, value) 
        {
            x@strand<-value 
            validObject(x) 
            return(x)
        })

setGeneric(name="position<-", def=function(object, value) standardGeneric("position<-"))
setReplaceMethod(f="position",
        signature="AlignedData", 
        definition=function(object, value) 
        {
            object@position<-value 
            validObject(object) 
            return(object)
        })

setGeneric(name="mapq<-", def=function(object, value) standardGeneric("mapq<-"))
setReplaceMethod(f="mapq",
        signature="AlignedData", 
        definition=function(object, value) 
        {
            object@mapq<-value 
            validObject(object) 
            return(object)
        })

setGeneric(name="isize<-", def=function(object, value) standardGeneric("isize<-"))
setReplaceMethod(f="isize",
        signature="AlignedData", 
        definition=function(object, value) 
        {
            object@isize<-value 
            validObject(object) 
            return(object)
        })

setGeneric(name="qwidth<-", def=function(x, value) standardGeneric("qwidth<-"))
setReplaceMethod(f="qwidth",
        signature="AlignedData", 
        definition=function(x, value) 
        {
            x@qwidth<-value
            validObject(x)
            return(x)
        })

setGeneric(name="flag<-", def=function(object, value) standardGeneric("flag<-"))
setReplaceMethod(f="flag",
        signature="AlignedData", 
        definition=function(object, value) 
        {object@flag<-value 
            validObject(object) 
            return(object)
        })

setGeneric(name="mrnm<-", def=function(object, value) standardGeneric("mrnm<-"))
setReplaceMethod(f="mrnm",
        signature="AlignedData", 
        definition=function(object, value) 
        {
            object@mrnm<-value 
            validObject(object) 
            return(object)
        })

setGeneric(name="notAligned<-", def=function(object, value) standardGeneric("notAligned<-"))
setReplaceMethod(f="notAligned",
        signature="AlignedData", 
        definition=function(object, value) 
        {
            object@notAligned<-value 
            validObject(object) 
            return(object)
        })

setGeneric(name="pairedEnds<-", def=function(object, value) standardGeneric("pairedEnds<-"))
setReplaceMethod(f="pairedEnds",
        signature="AlignedData", 
        definition=function(object, value) 
        {object@pairedEnds<-value 
            validObject(object) 
            return(object)
        })

setGeneric(name="weight<-", def=function(object, value) standardGeneric("weight<-"))
setReplaceMethod(f="weight",
        signature="AlignedData", 
        definition=function(object, value) 
        {
            object@weight<-value 
            validObject(object) 
            return(object)
        })



### Subsetters

#setGeneric(name="[", def=function(object, value) {standardGeneric("[")})
setMethod(f="[",
        signature="AlignedData", 
        definition=function(x,i,j,drop=FALSE) 
        {
            if(length(i)==0) return(initialize(x))
            
            testLength <- function(param, index) {return(if(length(param)>1) param[index] else param)} # This function is necessary not to try indexing vector of length 0 (optional slots, which returns NA...)
            return(initialize(x, 
                            readID=testLength(readID(x),i), 
                            seqnames=testLength(seqnames(x),i), # initialize eventually reduce levels 
                            strand=testLength(strand(x),i), # initialize eventually reduce levels
                            position=testLength(position(x),i), 
                            mapq=testLength(mapq(x),i), 
                            isize=testLength(isize(x),i), 
                            qwidth=testLength(qwidth(x),i), 
                            flag=testLength(flag(x),i), 
                            mrnm=testLength(mrnm(x),i), 
                            notAligned=notAligned(x), 
                            pairedEnds=pairedEnds(x),
                            weight=testLength(weight(x),i),
                            drop=drop))
        })



### Other member functions


setGeneric(name="dropUndefinedStrand", def=function(object, quiet=TRUE) standardGeneric("dropUndefinedStrand"))
setMethod(f="dropUndefinedStrand", signature="AlignedData", definition=function(object, quiet=TRUE)
        {
            undefinedStrand <- is.na(object@strand) | object@strand=="*"
            
            if((!quiet) && any(undefinedStrand)) 
            {
                cat("\n Filtering read(s) with undefined strand (", sum(undefinedStrand), ")", sep="")
                return(object[!undefinedStrand])
            }
            return(object)
        })

setGeneric(name="sortByPairs", def=function(object, quiet=TRUE) standardGeneric("sortByPairs"))
setMethod(f="sortByPairs", signature="AlignedData", definition=function(object, quiet=TRUE)
        {
            if(pairedEnds(object))
            {
                if(!quiet) cat("\n Sorting aligned data by pairs")
                indexesOrder <- order(readID(object), bitAnd(flag(object), 64)) # tells to put the 'first in pair' first, change to 128 for 'second in pair' first
                return(object[indexesOrder]) # be careful with eventual atomic values in the list that could be lost by subsetting
            }
            else
            {
                warning("Trying to sort by Pairs a dataset not declared as Paired-ends, ignoring...")
                return(object)
            }
        })

setGeneric(name="dropChromosomePattern", def=function(object, pattern, quiet=TRUE) standardGeneric("dropChromosomePattern"))
setMethod(f="dropChromosomePattern", signature="AlignedData", definition=function(object, pattern, quiet=TRUE)
        {
            if(!quiet) cat("\n Filtering reads matching chromosome pattern :", pattern)
            return(object[!grepl(pattern, seqnames(object))])
        })

setGeneric(name="filterInsertSize", def=function(object, rangeMin, rangeMax, includeLower=FALSE, quiet=TRUE) standardGeneric("filterInsertSize"))
setMethod(f="filterInsertSize", signature="AlignedData", definition=function(object, rangeMin, rangeMax, includeLower=FALSE, quiet=TRUE)
        {
            if(rangeMax<rangeMin) stop("Invalid rang of selection, rangeMax must be >= rangeMin")
            if(!quiet) cat("\n Filtering pairs by size range (min ", rangeMin, ", max ",rangeMax,")",sep="")
            if(length(isize(object))==0) 
            {
                stop("Trying to filter pairs by insert size but no insert size information is available")
            }
            
            if(includeLower)
            {
                return(object[(abs(isize(object))>=rangeMin) & (abs(isize(object))<=rangeMax)])
            }
            return(object[(abs(isize(object))>rangeMin) & (abs(isize(object))<=rangeMax)]) # operator '[' takes care of atomic values and reduces factors levels
        })

setGeneric(name="filterReadSize", def=function(object, rangeMin, rangeMax, includeLower=FALSE, quiet=TRUE) standardGeneric("filterReadSize"))
setMethod(f="filterReadSize", signature="AlignedData", definition=function(object, rangeMin, rangeMax, includeLower=FALSE, quiet=TRUE)
        {
            if(!quiet) cat("\n Filtering reads by size range (min ", rangeMin, ", max ",rangeMax,")",sep="")
            if(length(qwidth(object))==0)
            {
                stop("Trying to filter reads by read size but no read size information is available")
            }
            if(includeLower)
            {
                return(object[(qwidth(object)>=rangeMin) & (qwidth(object)<=rangeMax)])
            }
            return(object[(qwidth(object)>rangeMin) & (qwidth(object)<=rangeMax)]) # operator '[' takes care of atomic values and reduces factors levels
        })


# Generates a list of AlignedData objects, split by chromosomes, created as a list of environments
# Create environments in parent environment, to avoid them being erased at the end of the function, and return a list of them (~pointers)
#setGeneric(name="splitByChr_env", def=function(object, quiet=TRUE){standardGeneric("splitByChr_env")})
#setMethod(f="splitByChr_env", signature="AlignedData", definition=function(object, quiet=TRUE)
#        {
#            stop("Not implemented function")
#            
#            if(!quiet) cat("\n Splitting reads information by chromosomes")
#            splittedObject=split(object,chromosome(object))
#            
#            if(!quiet) cat("\n Encapsulating chromosomes in environments")
#            return(lapply(splittedObject,function(x)
#                    {
#                        capsule=new.env() # Create the new environment
#                        assign("value", x, envir=capsule) # Put the object in the environment
#                    }))
#            # Splitting the object in a list, each element being a new object containing information for a chromosome
#
#        })

# In case of paired-Ends, gives the indexes of orphans reads (the ones with no mate)
setGeneric(name="getOrphansIndexes", def=function(object, quiet=TRUE) standardGeneric("getOrphansIndexes"))
setMethod(f="getOrphansIndexes", signature="AlignedData", definition=function(object, quiet=TRUE)
        {
            if(pairedEnds(object))
            {
                if(!quiet) cat("\n Searching for orphans reads, (incomplete pairs)")
                return(!(readID(object) %in% readID(object)[duplicated(readID(object))]))
            }
            else
            {
                warning("Trying to retrieve orphans reads from a dataset not declared as Paired-ends, ignoring...")
                return(NULL)
            }
        })

# There is some room for interpretation on how the alignemnts must be reported, a lot of different situations can appear
# Hence some aligners (like tophat) can output some reads with complex properties (specially in complex situation such as paired ends)
# This function filters the reads that should not be processed by the pipeline (ie. 'mate unmapped', 'both reads on same strand')
# NOTE 1 : it does not rely on 'proper pair' flag which interpretation can vary among aligners (insert size for instance), it is up
# to the user to filter what is a proper pair either before using the pipeline or using the dedicated options of the pipeline (ignoreInsertsOver)
# NOTE 2 : the reads with 'mate unmapped' flag set are reported twice, both with the same coordinates
setGeneric(name="cleanNonSimplePairs", def=function(object, quiet=TRUE) standardGeneric("cleanNonSimplePairs"))
setMethod(f="cleanNonSimplePairs", signature="AlignedData", definition=function(object, quiet=TRUE)
        {
            if(pairedEnds(object))
            {
                if(!quiet) cat("\n Filtering pairs with  ...")
            }
            else
            {
                warning("cleanNonSimplePairs is inapropriate for a dataset not declared as Paired-ends, ignoring...")
                return(NULL)
            }
            
            #### Step 1. treat separately the "mate unmapped" reads which can be reported twice (for some of them both reads of the pair are reported with same coordinates and no isize, which can be problematic if one want to detect orphans based on duplicated readIDs)
            
            # get the read with "mate unmapped" flag set and remove the second match of each, they will eventually be considered as orphans in further steps of the pipeline
            pairsWithMateUnmappedFlag=( duplicated(readID(object)) & as.logical(bitAnd(flag(object), 8)) )
            if(sum(pairsWithMateUnmappedFlag)>0)
            {
                if(!quiet) cat("\n     Number of PAIRS reported with 'mate unmapped' flag set :",sum(pairsWithMateUnmappedFlag),"! Filtering (orphanize)...")
                object <- object[!pairsWithMateUnmappedFlag]
            }
            else
            {
                if(!quiet) cat("\n     No pairs reported with 'mate unmapped' flag set")
            }
            
            #### Step 2. "unlink" the pairs where both reads align on the same strand (very unlikely to be a valid pair (as opposed to what some of aligner report in the flag of some reads, they probably base their call on reads distance only)
            #### they will be considered as orphans un further steps of the pipeline (just have to assign a different readID for both reads of the pair)
            
            # get the pairs with both reads on the same strand (flag read reverse strand == flag mate reverse strand)
            readsSameStrand <- ( as.logical(bitAnd(flag(object), 16)) == as.logical(bitAnd(flag(object), 32)) )
            # get one of both reads in the pair
            readsToUnlink <- duplicated(readID(object)) & readsSameStrand
            if(sum(readsToUnlink)>0)
            {
                if(!quiet) cat("\n     Number of PAIRS reported with both reads on the same strand :",sum(readsToUnlink),"! Filtering (orphanize)...")
                # generate new readIDs for unlinking the pairs (as Pasha consider the pairs based on the readIDs duplication)
                maxreadID <- max(readID(object))
                newreadIDs <- maxreadID+(1:sum(readsToUnlink))
                # reaffect readIDs to unlink pairs (make them orphans)
                readID(object)[readsToUnlink] <- newreadIDs
            }
            else
            {
                if(!quiet) cat("\n     No pairs reported with both reads on the same strand")
            }
            
            return(object)
            
        })

# In case of paired-Ends, checks for reads consistency (all are correct pairs)
setGeneric(name="checkPairsOK", def=function(object) standardGeneric("checkPairsOK"))
setMethod(f="checkPairsOK", signature="AlignedData", definition=function(object)
        {
            if(pairedEnds(object))
            {
                cat("\n Checking pairs consistency")
                
                # Check that there is indeed information about pairs for all reads in the sam (using the first bit of the flag)
                if(length(object)%%2) 
                {
                    cat("\n     The number of reads in your file is not even... Very unlikely that the file have a mate for each read")
                    return(FALSE)
                }
                if(!all(as.logical(bitAnd(1,flag(object))))) 
                {
                    cat("\n     The experiment has been declared as paired end but at least some reads seem not to carry 'paired' information (check first bit of 'flag' in SAM file)")
                    return(FALSE)
                }
                
                cat("\n     OK")
                
                # Checks that all qnames (sequences ids) are periodic
                cat("\n Check that sorting seems ok...")
                even <- seq(2,length(object), by=2)
                odd <- seq(1,length(object), by=2)
                #if(!all( id(alignedDataList[["alignedObject"]][even]) == id(alignedDataList[["alignedObject"]][odd]) ))
                #Because its BStringSet and not characters we cannot use string comparison !
                if(!(all( readID(object)[even] %in% readID(object)[odd]) & all(readID(object)[odd] %in% readID(object)[even])))
                {
                    cat("\n     The reads and mates don't seem to come in pairs, please check that each read in your file has its mate in next/previous line or consider removing orphans and using the sort option")
                    return(FALSE)
                }
                
                cat("\n     OK")
                
                return(TRUE)
            }
            else
            {
                warning("Trying to check pairs consistency from a dataset not declared as Paired-ends, ignoring...")
                return(TRUE)
            }
        })

setGeneric(name="normalizeChrNames", def=function(object, chrPrefix, chrSuffix) standardGeneric("normalizeChrNames"))
setMethod(f="normalizeChrNames", signature="AlignedData", definition=function(object, chrPrefix="chr", chrSuffix="")
        {
            levels(seqnames(object))=paste("chr", gsub(paste(chrPrefix, chrSuffix,sep="|"),"", levels(seqnames(object))),sep="")
            return(object)
        })



### Public Constructor(s)

# Constructor from a file
readAlignedData = function(folderName, fileName, fileType="BAM", pairedEnds=FALSE)
{
    alignedDataObject <- NULL
    
    if((fileType!="BAM") & pairedEnds) stop("Paired-end reads support is limited to BAM format.")
    
    if(fileType=="BAM")
    {
        filePath <- if(nchar(folderName)>0) file.path(folderName, fileName) else fileName 
        if((!file.exists(filePath)) || (file.info(filePath)$isdir)) stop("The BAM filename cannot be found or is not a valid file...")
        
        # Use the BamFile 'connection' to handle the file reading (required for passing asmates parameter and ge the extra fields 'groupid' and 'mate_status')
        myFile <- BamFile(filePath, asMates=pairedEnds)
                
        # Define the flags to be respected by reads to be selected ('isNotPrimaryRead' has been deprecated at some point)
        bamFlags <- if(packageVersion("Rsamtools")>='1.20.4') scanBamFlag(isSecondaryAlignment=FALSE) else scanBamFlag(isNotPrimaryRead=FALSE)
        
        # Read the records, removing the non-simple CIGAR and the secondary alignments
        dataBAM <- scanBam(file=myFile, param=ScanBamParam(simpleCigar=TRUE, flag=, what=c("mapq", "pos", "strand", "rname", "qname", "isize", "qwidth", "flag", "mrnm", "groupid", "mate_status")))[[1]]
        
        # Finds the reads that were reported as aligned by the aligner (flag bit 0x4 == 0)
        isAligned <- !bitAnd(dataBAM[["flag"]], 4)
        
        # Removes the not aligned reads from the data
        dataBAM <- lapply(dataBAM, "[", isAligned)

        alignedDataObject <- new("AlignedData", readID=if(pairedEnds) dataBAM[["groupid"]] else as.integer(as.factor(dataBAM[["qname"]])), seqnames=dataBAM[["rname"]], strand=dataBAM[["strand"]], position=dataBAM[["pos"]], mapq=dataBAM[["mapq"]], isize=dataBAM[["isize"]], qwidth=dataBAM[["qwidth"]], flag=dataBAM[["flag"]], mrnm=dataBAM[["mrnm"]], notAligned=sum(!isAligned), pairedEnds=pairedEnds)
        
        
    }
    else if(fileType=="BED")
    {
        # TODO
        stop("Input file format not yet supported")
        #return(NULL)
    }
    else if(fileType=="GFF")
    {
        # TODO
        stop("Input file format not yet supported")
        #return(NULL)
    }
    else # Solexa and other types painless handled by readAligned from ShortRead
    {
        # Using the ShortRead parser to create a shortread object and extract its information to create a AlignedData object used in the program
        alignedData <- readAligned(dirPath=folderName, pattern=paste("^",fileName,"$",sep=""), type=fileType)
        
        alignedDataObject <- new("AlignedData", readID=1:length(alignedData), seqnames=as.character(seqnames(alignedData)), strand=as.character(strand(alignedData)), position=position(alignedData), mapq=slot(mapq(alignedData),"quality"), qwidth=width(sread(alignedData)), notAligned=as.integer(0), pairedEnds=FALSE)
    }
    
    return(alignedDataObject)
}



# Constructor from a file, reads output from multireads pipeline (reads that aligned on multiple positions)
.readMultipleAlignedData = function(fileName)
{
    multiLocData <- read.table(fileName, header=FALSE, quote="", comment.char="#", sep="\t")
    # Keless compatible format (no read length information)
    columnsNames <- c("chr", "strand", "pos", "weight")
    # Complete format
    if(ncol(multiLocData)==5)
    {
        columnsNames <- c("chr", "strand", "pos", "readLength", "weight")
    }
    colnames(multiLocData) <- columnsNames
    
    
    readLength <- integer(0)
    if(!is.null(multiLocData[["readLength"]])) 
    {
        readLength <- multiLocData[["readLength"]]
    }
    
    return(new("AlignedData", 
                    readID=1:nrow(multiLocData), 
                    seqnames=multiLocData[["chr"]], 
                    strand=multiLocData[["strand"]], 
                    position=multiLocData[["pos"]], 
                    qwidth=readLength, 
                    notAligned=as.integer(0), 
                    pairedEnds=FALSE, 
                    weight=multiLocData[["weight"]]))
    
}
