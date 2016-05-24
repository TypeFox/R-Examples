
library("sprint")
library("Biostrings")
library("BSgenome.Celegans.UCSC.ce2")
library("ff")
library("RUnit")
library("stringdist")
library("bigmemory")

## Randomly extract 10000 40-mers from C.elegans chrI:
extractRandomReads <- function(subject, nread, readlength)
{
	if (!is.integer(readlength))
	readlength <- as.integer(readlength)
	start <- sample(length(subject) - readlength + 1L, nread,
					replace=TRUE)
	DNAStringSet(subject, start=start, width=readlength)
}

args <- commandArgs(trailingOnly = TRUE)
nreads <- as.numeric(args[1])
print("number of reads")
print(nreads)

length <- as.numeric(args[2])
print("length")
print(length)

filename_ <- args[3]
print("filename")
print(filename_)

# Change this to test scaling.
rndreads <- extractRandomReads(Celegans$chrI, nreads, length)

stime_sprint <- proc.time()["elapsed"]
actual_result <- pstringDist(x=rndreads, method="h", filename=filename_)
etime_sprint <- proc.time()["elapsed"]

expected_result <- stringdistmatrix(rndreads, rndreads, method="h")


size <- nreads
desc5 <- new("big.matrix.descriptor"
, description = structure(list(sharedType = "FileBacked", filename = filename_, 
totalRows = size, totalCols = size, rowOffset = c(0, size), colOffset = c(0, 
size), nrow = size, ncol = size, rowNames = NULL, colNames = NULL, 
type = "integer", separated = FALSE), .Names = c("sharedType", 
"filename", "totalRows", "totalCols", "rowOffset", "colOffset", 
"nrow", "ncol", "rowNames", "colNames", "type", "separated"))
)

bigmemory_result <- attach.big.matrix(desc5)


checkTrue(all.equal(expected_result, actual_result[,], check.attributes=FALSE), "pstringDist and stringDist with list of strings.")
#checkEquals(dimnames(expected_result), dimnames(actual_result[,]), "Test labels on dist")
checkTrue(all.equal(expected_result, bigmemory_result[,], check.attributes=FALSE), "bigmemory check.")

print(paste("SPRINT pstringDist time: ")); print(paste(etime_sprint-stime_sprint))

print(system(paste("ls -lh|grep ",filename_)))
# system(paste("rm ",filename_))

