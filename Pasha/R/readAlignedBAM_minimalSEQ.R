# library(ShortRead)
# library(bitops)


# This function generates a readAligned (see ShortRead library) Object from a BAM file.
# However it ignores the sequences and quality strings and replaces it by a single artificial nucleotide.
# This method helps to overcome the 2^31 limit of nucleotides when reading big BAM files if the final user
# doesn't care about the sequence content. This could also help to limit the memory footprint of the object.
.readAlignedBAM_minimalSEQ <- function(fileName)
{
    if(length(fileName)!=1) stop("fileName argument length is not 1")
    if((!file.exists(fileName)) || (file.info(fileName)$isdir)) stop("The filename is not a valid file")
    
    # Read the BAM file, here the outer resulting list has length one since no 'which' is defined in scanBamParam
    # TODO : handle 'which' argument
    dataBAM <- scanBam(file=fileName, param=ScanBamParam(what=c("mapq", "pos", "strand", "rname", "qname", "isize", "qwidth", "flag", "mrnm")))[[1]]
    
    # Finds the reads that were reported as aligned by the aligner (flag bit 0x4 == 0)
    isAligned <- !bitAnd(dataBAM[["flag"]], 4)
    nbAligned <- sum(isAligned)
    
    # Removes the not aligned reads from the data
    dataBAM <- lapply(dataBAM, "[", isAligned)
    
    # Builds 'readAligned' object with aligned reads, filling only partial information (no sequence nor quality and arbitrary ID)
    alignedObject  <-  ShortRead::AlignedRead(sread=Biostrings::DNAStringSet(rep("N", nbAligned)),
            id=Biostrings::BStringSet(dataBAM[["qname"]]),
            quality=ShortRead::FastqQuality(rep("%", nbAligned)),
            seqnames=dataBAM[["rname"]],
            strand=GenomicRanges::strand(dataBAM[["strand"]]),
            position=as.integer(dataBAM[["pos"]]),
            mapq=ShortRead::NumericQuality(dataBAM[["mapq"]])
    )
    
    return(list(alignedObject=alignedObject, isize=dataBAM[["isize"]], qwidth=dataBAM[["qwidth"]], flag=dataBAM[["flag"]], mrnm=dataBAM[["mrnm"]], notAligned=sum(!isAligned)))
}


.readAlignedBAM_minimalSEQ_NOSHORTREAD <- function(fileName, simpleCigar=TRUE)
{
    if(length(fileName)!=1) stop("fileName argument length is not 1")
    if((!file.exists(fileName)) || (file.info(fileName)$isdir)) stop("The filename is not a valid file")
    
    # Read the BAM file, here the outer resulting list has length one since no 'which' is defined in scanBamParam
    # TODO : handle 'which' argument
    dataBAM <- scanBam(file=fileName, param=ScanBamParam(simpleCigar=simpleCigar, what=c("mapq", "pos", "strand", "rname", "qname", "isize", "qwidth", "flag", "mrnm")))[[1]]
    
    # Finds the reads that were reported as aligned by the aligner (flag bit 0x4 == 0)
    isAligned <- !bitAnd(dataBAM[["flag"]], 4)
    
    # Removes the not aligned reads from the data
    dataBAM <- lapply(dataBAM, "[", isAligned)
    
    return(list(readID=dataBAM[["qname"]], seqnames=dataBAM[["rname"]], strand=strand(dataBAM[["strand"]]), position=as.integer(dataBAM[["pos"]]), mapq=dataBAM[["mapq"]], isize=dataBAM[["isize"]], qwidth=dataBAM[["qwidth"]], flag=dataBAM[["flag"]], mrnm=dataBAM[["mrnm"]], notAligned=sum(!isAligned)))
}
