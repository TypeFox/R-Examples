### R code from vignette source 'Pasha.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: testAlignedDataClass
###################################################
library(Pasha);

# Get the location of the test file bundled in the package    
testFile <- system.file("extdata", "embedDataTest.bam",package="Pasha");

# Load the file in an object
myData <- readAlignedData(
folderName <- dirname(testFile),
fileName <- basename(testFile),
fileType <- "BAM",
pairedEnds <- TRUE);

cat("\nNumber of reads in the loaded file :", length(myData), "\n");

# Remove reads that don't have a mate (orphans)
myData <- myData[!getOrphansIndexes(myData, quiet=FALSE)];

# Checking that paired-end representation is compatible with internal
# class representation
if(!checkPairsOK(myData))
{
    myData <- sortByPairs(myData, quiet=FALSE);
}

# Split the dataset in a list by chromosomes
myData_ChrList <- split(myData, seqnames(myData));

# Plot the number of reads per chromosome in the dataset
barplot(sapply(myData_ChrList, length),
        main="Number of reads by chromosome",
        col=terrain.colors(length(myData_ChrList)));


