readBed <-
function(bed.file=NULL) {
    if(is.null(bed.file)) bed.file <- file.choose()
    cat("Reading", bed.file , "...\n")
    bed.header <- readLines(bed.file, n=10)
    bed.start <- min(grep("chr", bed.header))
    bed.table <- read.table(bed.file, skip=bed.start-1, header=FALSE, sep="\t")
    colnames(bed.table)[1:3] <- c("chrom", "chromStart", "chromEnd")
    return(bed.table)
}
