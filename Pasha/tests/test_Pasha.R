
library(Pasha)

###########################################
## 1. Classic tests for the class
###########################################

# Get the BAM file used in RSamTools as example
BAMfile <- system.file("extdata","ex1.bam",package="Rsamtools")
# Testing the public constructor of AlignedData class (from Pasha)
alignedDataObject <- readAlignedData(folderName="" ,
                                     fileName=BAMfile, 
                                     fileType="BAM", 
                                     pairedEnds=T)
print(alignedDataObject)
                             
# Rename the chromosomes from seqX to chrX
alignedDataObject <- normalizeChrNames(alignedDataObject,
                                       chrPrefix="seq",
                                       chrSuffix="")
print(alignedDataObject)
                               
                               
#### Splitting object (indirect use of '[') 
# Split the object and the associated data by chromosomes 
splitListalignedDataObject <- split(alignedDataObject, 
                                    seqnames(alignedDataObject))
# check that the object has correctly been split
print(names(splitListalignedDataObject))
print(sapply(lapply(splitListalignedDataObject, seqnames), unique))


#### Orphans and pairs consistency
print(length(alignedDataObject))
# Identify the reads with no mate and remove them
orphansReadsIndex <- getOrphansIndexes(alignedDataObject)
print(paste("Nb of orphan reads :", sum(orphansReadsIndex)))
alignedDataObject <- alignedDataObject[!orphansReadsIndex]
print(length(alignedDataObject))

# Check that pairs are correctly represented in the object
# ie. mate pairs are consecutive
print(checkPairsOK(alignedDataObject))
# Sort by pairs
alignedDataObject <- sortByPairs(alignedDataObject, quiet=FALSE)
# Check pairs consistency again
print(checkPairsOK(alignedDataObject))


#### Filtering

## Chromosome patterns
# (filter reads on chromosomes containing the character '1')
print(alignedDataObject)
resultFiltering <- dropChromosomePattern(alignedDataObject, 
                                           pattern="1")
print(resultFiltering)

## Range of insert size
print(table(isize(alignedDataObject)))
resultFiltering <- filterInsertSize(alignedDataObject, 
                                    rangeMin=190, 
                                    rangeMax=195)
print(table(isize(resultFiltering)))
print(resultFiltering)

## Reads sequence size
print(table(qwidth(alignedDataObject)))
resultFiltering <- filterReadSize(alignedDataObject, 
                                  rangeMin=35, 
                                  rangeMax=50)
print(table(qwidth(resultFiltering)))
print(resultFiltering)
                          
## Reads with atypic characteristics for chromatin alignments
# Mate unmapped flag set, or mates aligned on same strand
resultFiltering <- cleanNonSimplePairs(alignedDataObject)
# Reads with a strand value different from expected '+' or '-'
resultFiltering <- dropUndefinedStrand(resultFiltering, quiet=FALSE)




###########################################
## 2. Launch Pasha package functional testing (regular test, a complete one can be started later manually by user, but is time consuming)
###########################################

Pasha:::.testFunctional(tempdir(), testType="regular", verbose=FALSE)

###########################################
## 3. Launch Pasha package functional multiread testing 
###########################################

Pasha:::.testFunctionalMultiread(tempdir(), verbose=FALSE)
