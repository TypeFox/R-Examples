
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#   File   : rbamtools.r
#   Date   : 12.Mar.2012
#   Sam    : Samtools downloaded September 7, 2011. Format: v1.4-r985
#   Content: R-Source for package rbamtools
#   Version: 2.10.8
#   Author : W. Kaisers
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#  Changelog
#  30.Okt.12  [initialize.gapList] Made printout message optional (verbose)
#  31.Okt.12  [bamRange] Included test for initialized index
#  31.Okt.12  Check for open reader in getHeader, getHeaderText, getRefCount
#  01.Nov.12  [get_const_next_align] added to correct memory leak.
#  08.Nov.12  Reading and writing big bamRanges (pure C, no R) valgrind
#                                                                checked.
#  09.Nov.12  [bamCopy.bamReader] Added which allows refwise copying.
#  31.Dec.12  gapSiteList class added
#  11.Jan.13  bamGapList class added
#  06.Feb.13  First successful test of bamGapList on 36 BAM-files 
#                                                       (871.926/sec)
#  20.Feb.13  Fixed Error in merge.bamGapList
#  27.Feb.13  Renamed createIndex -> create.index and loadIndex -> load.index
#             and bamSiteList -> bamGapList
#  18.Apr.13  Corrected some memory leaks in C-Code as reported by Brian Ripley
#  22.Apr.13  Added (read-) name and revstrand to as.data.frame.bamRange
#                                               (as proposed by Ander Muniategui)
#  11.Jun.13  Added reader2fastq and range2fastq functions
#                                               (2.5.3, valgrind tested)
#  11.Jun.13  Changed signature for bamSave: added refid argument
#               (needed to prevent samtools crashes when creating BAM files
#               with single align regions and appropriate refSeqDict entries)
#             (2.5.4, valgrind tested)
#  12.Jun.13  Added extractRanges function (2.5.5)
#  21.Jun.13  Added bamAlign function (2.5.6)
#  01.Jul.13  Added bamCount function (2.5.8)
#  02.Jul.13  Added bamCountAll function, valgrind tested (2.5.9)
#  18.Jul.13  Changed 'nGapAligns' to 'nAlignGaps (2.5.10)
#                                               nGapAligns deprecated!
#  24.JUl.13  Added alignQual function, valgrind tested (2.5.11)
#  28.Jul.13  Added alignDepth function, valgrind tested (2.5.12)
#  13.Aug.13  Added countTextLines function, valgrind tested (2.6.1)
#  26.Aug.13  Removed "coerce" from Namespace declaration
#  02.Sep.13  Changed "cat" to "message"
#  10.Jun.14  On CRAN after correction of "Mis-alignment errors"
#  19.Jun.14  Re-introduction of changes after resetting to 2.7.0
#                                                       due to internal errors.
#  21.Jul.14  Updated plotAlignDepth 
#  14.Jul.14  Added test directory
#  28.Jul.14  Added NEWS and ChangeLog file
#  29.Sep.14  Added support for DS segment in headerProgram (@PG)
#               Added support for Supplementary alignmnet FLAG
#               Corrected error in resetting FLAG values
#               Replaced rand() by runif() in ksort.h
#               Enclosed reader2fastq example in \dontrun{}
#  03.Nov.14  Changed nAligns data type to unsigned long long int
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

.onUnload <- function(libpath) { library.dynam.unload("rbamtools", libpath) }



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#   Declaration of generics
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

# Reader associated Generics

setGeneric("filename", function(object) standardGeneric("filename"))

setGeneric("isOpen", function(con, rw="") standardGeneric("isOpen"))

setGeneric("bamClose", function(object) standardGeneric("bamClose"))

setGeneric("rewind", function(object) standardGeneric("rewind"))

setGeneric("getHeader",function(object) standardGeneric("getHeader"))

setGeneric("getRefCount",function(object) standardGeneric("getRefCount"))

setGeneric("getRefData",function(object) standardGeneric("getRefData"))

setGeneric("getRefCoords",function(object, sn) standardGeneric("getRefCoords"))

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#  Replacement for (deprecated) functions -> .Defunct
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
setGeneric("createIndex",function(object,idx_filename)
    standardGeneric("createIndex"))

setGeneric("loadIndex", function(object, filename)
    standardGeneric("loadIndex"))

setGeneric("indexInitialized", function(object) 
    standardGeneric("indexInitialized"))

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#  Soon deprecated functions (for consistency reasons)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

setGeneric("create.index",function(object, idx_filename)
    standardGeneric("create.index"))

setGeneric("load.index",function(object, filename)
    standardGeneric("load.index"))

setGeneric("index.initialized",function(object) 
    standardGeneric("index.initialized"))

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

setGeneric("bamSort",function(object, prefix="sorted",
                              byName=FALSE, maxmem=1e+9,
                              path=dirname(filename(object))) standardGeneric("bamSort"))

setGeneric("reader2fastq",function(object, filename, which, append=FALSE)
    standardGeneric("reader2fastq"))

setGeneric("readerToFastq", function(object, filename, which, append=FALSE)
    standardGeneric("readerToFastq"))

setGeneric("bamCopy", function(object, writer, refids, verbose=FALSE)
    standardGeneric("bamCopy"))

setGeneric("extractRanges",function(object, ranges, filename, complex=FALSE,
                                    header, idxname) standardGeneric("extractRanges"))

setGeneric("bamCount", function(object, coords) standardGeneric("bamCount"))

setGeneric("bamCountAll", function(object, verbose=FALSE)
    standardGeneric("bamCountAll"))

setGeneric("nucStats", function(object, ...) standardGeneric("nucStats"))


# generic for bamReader and bamRange
setGeneric("getNextAlign", function(object) standardGeneric("getNextAlign"))

# generic for bamWriter and bamRange
setGeneric("bamSave", function(object, ...) standardGeneric("bamSave"))

# Generic for conversion into list
setGeneric("as.list", function(x, ...) standardGeneric("as.list"))

# Generic for retrieving RefData string from Objects
setGeneric("getHeaderText", function(object, delim="\n")
    standardGeneric("getHeaderText"))

# Generic for Reading member from object list
setGeneric("getVal", function(object, member) standardGeneric("getVal"))

# Generic for Writing member to object list
setGeneric("setVal", function(object, members, values)
    standardGeneric("setVal"))

# Generic fora adding read group to header
setGeneric("addReadGroup",function(object, l) standardGeneric("addReadGroup"))

# Generic for retrieving of list size
setGeneric("size", function(object) standardGeneric("size"))

# Generic for retrieving Nr of aligns in BAM region from gapList
setGeneric("nAligns", function(object) standardGeneric("nAligns"))

# Generic for retrieving Nr of gapped-aligns in BAM region from gapList
setGeneric("nAlignGaps", function(object) standardGeneric("nAlignGaps"))

# Generic for reading gapLists (align gaps) from bamReader
setGeneric("gapList", function(object, coords) standardGeneric("gapList"))

# Generic for reading gapSiteList (merged align gap sites) from bamReader
setGeneric("siteList", function(object, coords) standardGeneric("siteList"))

# Generic for reading bamGapList (merged align gap sites for whole bam-files)
# from bamReader
setGeneric("bamGapList", function(object) standardGeneric("bamGapList"))

# Generic for retrieving quality values
setGeneric("getQualDf", function(object, prob=FALSE, ...)
    standardGeneric("getQualDf"))

# Generic for retrieving quantile values from (phred) 
# quality tables (used for plotQualQuant)
setGeneric("getQualQuantiles", function(object, quantiles, ...)
    standardGeneric("getQualQuantiles"))

# Generic for plotting of (phred) quality quantiles.
setGeneric("plotQualQuant", function(object)
    standardGeneric("plotQualQuant"))

# bamHeader related generics
setGeneric("bamWriter", function(x, filename) standardGeneric("bamWriter"))

# refSeqDict related generics
setGeneric("removeSeqs", function(x, rows) standardGeneric("removeSeqs"))

setGeneric("addSeq",
           function(object, SN, LN, AS="", M5=0, SP="", UR="")
               standardGeneric("addSeq"))

setGeneric("head", function(x, ...) standardGeneric("head"))
setGeneric("tail", function(x, ...) standardGeneric("tail"))

# bamHeaderText related generics
setGeneric("headerLine", function(object) standardGeneric("headerLine"))

setGeneric("refSeqDict", function(object) standardGeneric("refSeqDict"))

setGeneric("headerReadGroup", function(object)
    standardGeneric("headerReadGroup"))

setGeneric("headerProgram", function(object)
    standardGeneric("headerProgram"))

setGeneric("headerLine<-", function(object,value)
    standardGeneric("headerLine<-"))

setGeneric("refSeqDict<-", function(object, value)
    standardGeneric("refSeqDict<-"))

setGeneric("headerReadGroup<-", function(object, value)
    standardGeneric("headerReadGroup<-"))

setGeneric("headerProgram<-", function(object, value)
    standardGeneric("headerProgram<-"))

setGeneric("bamHeader", function(object)
    standardGeneric("bamHeader"))

# gapSiteList related generics
setGeneric("refID", function(object) standardGeneric("refID"))

# bamRange related generics
setGeneric("getCoords", function(object) standardGeneric("getCoords"))

setGeneric("getParams", function(object) standardGeneric("getParams"))

setGeneric("getSeqLen", function(object) standardGeneric("getSeqLen"))

setGeneric("getRefName", function(object) standardGeneric("getRefName"))

setGeneric("getAlignRange", function(object) standardGeneric("getAlignRange"))

setGeneric("getPrevAlign", function(object) standardGeneric("getPrevAlign"))

setGeneric("stepNextAlign", function(object) standardGeneric("stepNextAlign"))

setGeneric("stepPrevAlign", function(object) standardGeneric("stepPrevAlign"))

setGeneric("push_back", function(object, value) standardGeneric("push_back"))

setGeneric("pop_back", function(object) standardGeneric("pop_back"))

setGeneric("push_front", function(object, value) standardGeneric("push_front"))

setGeneric("pop_front", function(object) standardGeneric("pop_front"))

setGeneric("writeCurrentAlign", 
           function(object, value) standardGeneric("writeCurrentAlign"))

setGeneric("insertPastCurrent", 
           function(object, value) standardGeneric("insertPastCurrent"))

setGeneric("insertPreCurrent", 
           function(object, value) standardGeneric("insertPreCurrent"))

setGeneric("moveCurrentAlign",
           function(object, target) standardGeneric("moveCurrentAlign"))

# Deprecated:
setGeneric("range2fastq", 
           function(object, filename, which, append=FALSE)
               standardGeneric("range2fastq"))

setGeneric("rangeToFastq", 
           function(object, filename, which, append=FALSE)
               standardGeneric("rangeToFastq"))

setGeneric("countNucs", function(object) standardGeneric("countNucs"))


# alignDepth related generics
setGeneric("alignDepth", 
           function(object, gap=FALSE) standardGeneric("alignDepth"))

setGeneric("getDepth", 
           function(object, named=FALSE) standardGeneric("getDepth"))

setGeneric("getPos", function(object) standardGeneric("getPos"))

setGeneric("plotAlignDepth",
           function(object, start=NULL, end=NULL, xlim=NULL,
                    main="Align Depth", xlab="Position", 
                    ylab="Align Depth",  transcript="",
                    strand=NULL , log="y", cex.main=2,
                    col="grey50", fill="grey90", grid=TRUE, 
                    box.col="grey20", box.border="grey80", ... )
               standardGeneric("plotAlignDepth")
)

# bamAlign related generics
setGeneric("name", function(object)
    standardGeneric("name"))

setGeneric("position", function(object)
    standardGeneric("position"))

setGeneric("nCigar", function(object)
    standardGeneric("nCigar"))

setGeneric("cigarData", function(object)
    standardGeneric("cigarData"))

setGeneric("mateRefID", function(object)
    standardGeneric("mateRefID"))

setGeneric("matePosition", function(object)
    standardGeneric("matePosition"))

setGeneric("insertSize", function(object)
    standardGeneric("insertSize"))

setGeneric("mapQuality", function(object)
    standardGeneric("mapQuality"))

setGeneric("alignSeq", function(object)
    standardGeneric("alignSeq"))

setGeneric("alignQual", function(object)
    standardGeneric("alignQual"))

setGeneric("alignQualVal", function(object)
    standardGeneric("alignQualVal"))

setGeneric("pcrORopt_duplicate", function(object)
    standardGeneric("pcrORopt_duplicate"))

setGeneric("pcrORopt_duplicate<-", function(object, value)
    standardGeneric("pcrORopt_duplicate<-"))

setGeneric("failedQC", function(object)
    standardGeneric("failedQC"))

setGeneric("failedQC<-", function(object, value)
    standardGeneric("failedQC<-"))

setGeneric("firstInPair", function(object) 
    standardGeneric("firstInPair"))

setGeneric("firstInPair<-", function(object ,value)
    standardGeneric("firstInPair<-"))

setGeneric("secondInPair", function(object) 
    standardGeneric("secondInPair"))

setGeneric("secondInPair<-", function(object, value)
    standardGeneric("secondInPair<-"))

setGeneric("unmapped", function(object)
    standardGeneric("unmapped"))

setGeneric("unmapped<-", function(object, value)
    standardGeneric("unmapped<-"))

setGeneric("mateUnmapped", function(object)
    standardGeneric("mateUnmapped"))

setGeneric("mateUnmapped<-", function(object, value)
    standardGeneric("mateUnmapped<-"))

setGeneric("reverseStrand", function(object)
    standardGeneric("reverseStrand"))

setGeneric("reverseStrand<-", function(object, value)
    standardGeneric("reverseStrand<-"))

setGeneric("mateReverseStrand", function(object)
    standardGeneric("mateReverseStrand"))

setGeneric("mateReverseStrand<-", function(object, value)
    standardGeneric("mateReverseStrand<-"))

setGeneric("paired", function(object) 
    standardGeneric("paired"))

setGeneric("paired<-", function(object, value)
    standardGeneric("paired<-"))

setGeneric("properPair", function(object)
    standardGeneric("properPair"))

setGeneric("properPair<-", function(object, value)
    standardGeneric("properPair<-"))

setGeneric("secondaryAlign", function(object)
    standardGeneric("secondaryAlign"))

setGeneric("secondaryAlign<-", function(object, value)
    standardGeneric("secondaryAlign<-"))

setGeneric("suppAlign", function(object)
    standardGeneric("suppAlign"))

setGeneric("suppAlign<-", function(object, value)
    standardGeneric("suppAlign<-"))

setGeneric("flag", function(object) 
    standardGeneric("flag"))

setGeneric("flag<-", function(object, value)
    standardGeneric("flag<-"))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# GenomePartition class
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

setGeneric("getSeqNr", function(object)
                                    standardGeneric("getSeqNr"))

setGeneric("genomePartition", function(object, genome, dist=10000) 
                                    standardGeneric("genomePartition"))

setGeneric("getGridAlignCounts", function(object) 
                                    standardGeneric("getGridAlignCounts"))

setGeneric("getAlignCounts", function(object) 
                                    standardGeneric("getAlignCounts"))

# src = source
setGeneric("countPartition", function(partition, src) 
                                    standardGeneric("countPartition"))

setGeneric("getFileTable", function(object) 
                                    standardGeneric("getFileTable"))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# End setGeneric
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#  static functions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

gc_content <- function(An, Cn, Gn, Tn)
{
    denom <- sum(An+Cn+Gn+Tn)
    if(denom==0)
        return(0)
    return(sum(Gn + Cn) / denom)
}
at_gc_ratio <- function(An, Cn, Gn, Tn)
{
    denom <- sum(Gn+Cn)
    if(denom==0)
        return(NA)
    return(sum(An + Tn) / denom)
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#  Declaration of classes
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#  File interacting classes:
#  bamReader, bamWriter
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

setClass("bamReader", representation(
    filename="character", reader="externalptr", index="externalptr",
    startpos="numeric"),
    validity=function(object)
    {return(ifelse(is.null(object@reader), FALSE, TRUE))})


setClass("bamWriter",
            representation(filename="character", writer="externalptr"),
            validity=function(object)
            {
                return(ifelse(is.null(object@writer), FALSE, TRUE))
            })

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#  Header section related classes:
#  bamHeader, bamHeaderText,
#  headerLine,
#  
#  refSectDict, headerReadGroup,
#  headerProgram
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

setClass("bamHeader",
        representation(header="externalptr"), validity=function(object)
{
    return(ifelse(is.null(object@header), FALSE, TRUE))
})

setClass("headerLine",
    representation(VN="character", SO="character"),
    validity=function(object)
    {
        if( (length(VN) == 1) & (length(SO) == 1) )
            return(TRUE)
        else
            return(FALSE)
    })


setClass("refSeqDict",
    representation(SN="character", LN="numeric",   AS="character",
                   M5="numeric",   SP="character", UR="character"))

setClass("headerReadGroup", representation(
            nrg="integer",          #  Number of read groups
            ID="character",         #  Read group identifier
            CN="character",         #  Name of sequencing center
            DS="character",         #  Description
            DT="character",         #  Date
            FO="character",         #  Flow order
            KS="character",         #  Array of nucleotide bases
            LB="character",         #  Library
            PG="character",         #  Programes used for processing
            PI="character",         #  Predicted median insert size
            PL="character",         #  Platform (ILLUMINA,...)
            PU="character",         #  Platform unit (e.g. lane code)
            SM="character",         #  Sample (Pool name)
            ntl="integer"           #  number of taglabs (= 12 (static))
            ),
                                validity=function(object) {return(TRUE)})

tagLabs <- c("ID", "CN", "DS", "DT", "FO", "KS", "LB", "PG",
             "PI", "PL", "PU", "SM")

setClass("headerProgram",
    representation(l="list"),
    validity=function(object) {return(TRUE)})

setClass("bamHeaderText",
        representation(head="headerLine", dict="refSeqDict",
                    group="headerReadGroup", prog="headerProgram", 
                    com="character"))

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#  Align related classes:
#  bamAlign, bamRange, alignDepth
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

setClass("bamAlign", 
    representation(align="externalptr"), validity=function(object)
{
    return(ifelse(is.null(object@align, FALSE, TRUE)))
})


setClass("bamRange", 
    representation(range="externalptr"), validity=function(object)
{
    return(ifelse(is.null(object@range), FALSE, TRUE))
})


setClass("alignDepth",
            representation(depth="integer", pos="integer",
            params="numeric", refname="character"))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#  Segment Count
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Sequence segments: Tiling a single reference sequence (chromosome)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

setClass("seqSegments",
        representation(
            # Sequence header
            seqid = "numeric",
            seqname = "character",
            seqlen = "numeric",
            # Position data
            margins = "numeric"
        ),
        validity = function(object) {
            if(length(object@seqid)!=1 || 
                length(object@seqname)!=1 ||
                length(object@seqlen) !=1)
            { 
                message("[validity.seqSegments] Sequence values must have length 1!")
                return(FALSE)
            }
            
            if(any(object@margins < 0))
            {
                message("[validity.seqSegments] No negative margin values allowed!")
                return(FALSE)
            }
            
            if(any(object@margins > object@seqlen))
            {
                message("[validity.seqSegments] Out of bound margins (> seqlen)!")
                return(FALSE)
            }
            
            return(TRUE)
        }
)

seqSegments <- function(bam=NULL, i=1)
{
    if(is.null(bam))
        return( new("seqSegments"))
    
    if(!is.character(bam))
        stop("bam must be character!")
    
    if(!file.exists(bam))
        stop("bam does not exist!")
    
    reader <- bamReader(bam)
    rd <- getRefData(reader)
    
    res <- new("seqSegments")
    res@seqid <- rd$ID[i]
    res@seqname <- rd$SN[i]
    res@seqlen <- rd$LN[i]
    
    res@margins <- c(1, rd$LN[i])
    return(res)
}

setGeneric("margins", function(object) standardGeneric("margins"))
setMethod("margins", "seqSegments", function(object){
    return(object@margins)
})

setGeneric("margins<-", function(object, value) standardGeneric("margins<-"))
setReplaceMethod("margins", "seqSegments", function(object, value){
    if(!is.numeric(value))
        stop("value must be numeric!")

    object@margins <- sort(unique(value))
    if(object@margins[1] < 0)
        stop("Negative margin values are not allowed!")
    
    if(object@margins[length(object@margins)] > object@seqlen[1])
        stop("Margin values out of bounds (> seqlen)!")
    
    return(object)
})

setGeneric("addMargins", function(object, value=NULL) standardGeneric("addMargins"))
setMethod("addMargins", "seqSegments", function(object, value=NULL){
    
    if(is.null(value))
        return(object)
    
    if(!is.numeric(value))
        stop("value must be numeric!")
    
    object@margins <- sort(unique(c(object@margins, value)))
    
    if(!validObject(object))
        stop("")

    return(object)
})


setMethod("show", "seqSegments", function(object)
{
    bm<-Sys.localeconv()[7]
    r <-"right"
    w <- 14
        
    cat("An object of class '", class(object), "'.\n", sep="")
    cat("Reference sequence:\n")
    cat("Seqid                 : ", format(object@seqid[1],   w=w, j=r), "\n")
    cat("Seqname               : ", format(object@seqname[1], w=w, j=r), "\n")
    cat("Seqlen                : ", format(
            format(object@seqlen[1], big.m=bm), w=w, j=r), "\n")
    
    nm <- length(object@margins)
    
    if(nm == 0)
    {
        cat("Segments number       : 0\n")
        return(invisible())
    }
    
    cat("Segments number       : ", format(
        format(nm - 1, big.m=bm),               w=w, j=r), "\n")
    
    cat("Segments left  margin : ", format(
        format(object@margins[1], big.m=bm),    w=w, j=r), "\n")
    
    cat("Segments right margin : ", format(
        format(object@margins[nm], big.m=bm),   w=w, j=r), "\n")
    
    # Set output width so that all numbers can be printed at equal width
    w <- ceiling(log10(object@margins[nm])) + 3
    if(nm > 1)
    {
        cat("Segments :\n")
        nd <- min((nm - 1), 6)
        for(i in 1:nd)
        {
            cat(i, ": [", 
                format(format(object@margins[i], big.m=bm), w=w, j=r),
                ", ",
                format(format(object@margins[i + 1], big.m=bm), w=w, j=r),
                ")\n",
                sep="")
        }
    }
    return(invisible())
})



setClass("rangeSegCount",
        representation(
            position = "integer",
            count = "integer",
            refname = "character",
            LN = "integer",
            coords = "numeric",
            complex = "logical"
        ),
        validity = function(object) { return(length(position)==length(count)) }
)


setMethod("show", "rangeSegCount", function(object)
{
    o <- object
    op <- o@coords
    n <- 6
    cat("An object of class '", class(o), "'.\n",sep="")
    bm<-Sys.localeconv()[7]
    w<-15
    r<-"right"
    cat("Refname : ", format(o@refname              , w=w, j=r), "\n", sep="")
    cat("Seqid   : ", format(format(op[1], big.m=bm), w=w, j=r), "\n", sep="")
    cat("LN      : ", format(format(o@LN,  big.m=bm), w=w, j=r), "\n", sep="")
    cat("qrBegin : ", format(format(op[2], big.m=bm), w=w, j=r), "\n", sep="")
    cat("qrEnd   : ", format(format(op[3], big.m=bm), w=w, j=r), "\n", sep="")
    cat("Complex : ", format(o@complex              , w=w, j=r), "\n", sep="")
    cat("Size    : ", format(format(length(o@position),  big.m=bm), w=w, j=r), "\n", sep="")
    cat("\n")
    print(data.frame(position=o@position[1:n], count=o@count[1:n]))
    return(invisible())
})



setMethod("size", "rangeSegCount", function(object) {return(length(object@position))})

as.data.frame.rangeSegCount <- function(x, row.names=NULL, optional=FALSE, ...)
{
    return(data.frame(position=x@position, 
                        count=x@count,
                        row.names=row.names))
}

setGeneric("rangeSegCount", function(object, 
                coords=NULL, segments=NULL, 
                complex=FALSE) standardGeneric("rangeSegCount"))

setMethod("rangeSegCount", "bamReader", 
    function(object, coords=NULL, segments=NULL, complex=FALSE)
    {
        if(!indexInitialized(object))
            stop("Reader must have initialized index! Use 'load.index'!")
        
        if(is.null(coords))
            stop("coords argument is not optional!")
        
        if(!is.numeric(coords))
            stop("coords must be numeric!")
        
        if(length(coords) != 3)
            stop("coords must have length 3")
        
        if(any(coords < 0))
            stop("coords must not be negative")
        
        coords <- as.numeric(coords)
        
        if(is.null(segments))
            stop("segments argument is not optional")
        
        if(!is.numeric(segments))
            stop("segments must be numeric")
        
        if(any(segments < 0))
            stop("segment values must not be negative")
        
        segments <- sort(as.integer(segments))
        
        if(!is.logical(complex))
            stop("complex must be logical")
        
        res <- .Call("bam_count_segment_aligns",
                     object@reader, object@index,
                     coords, segments, 
                     complex[1], PACKAGE="rbamtools")
        
        return(res)
        
    }
)


setGeneric("meltDownSegments", 
           function(object, factor=1) standardGeneric("meltDownSegments"))

setMethod("meltDownSegments", "rangeSegCount", function(object, factor=1)
{
    if(!is.numeric(factor))
        stop("factor must be numeric")
    
    if(length(factor) != 1)
        stop("factor must contain exactly one value")
    
    if(factor[1] < 1)
        stop("factor must be >= 1")
    
    factor <- as.integer(factor)
    
    res <- .Call("bam_count_segment_melt_down", 
                 object, factor[1], PACKAGE="rbamtools")
    return(res)
})


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#  Gap sites related classes:
#  gapList, gapSiteList
#  bamGapList
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

setClass("gapList", 
    representation(list="externalptr"), validity=function(object)
{
    return(ifelse(is.null(object@list), FALSE, TRUE))
})

setClass("gapSiteList",
    representation(list="externalptr"), validity=function(object)
{
    return(ifelse(is.null(object@list), FALSE, TRUE))
})

setClass("bamGapList",
        representation(list="externalptr", refdata="data.frame"),
        validity=function(object)
{ 
    return(ifelse(is.null(object@list), FALSE, TRUE))
})




# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#
#  bamReader
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#   Opening and closing a BAM-File for reading
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

setMethod(f="initialize", signature="bamReader",
            definition=function(.Object, filename)
{
        .Object@filename <- filename
        .Object@reader <- .Call("bam_reader_open",
                                    path.expand(filename), PACKAGE="rbamtools")
        
        .Object@startpos <- .Call("bam_reader_tell", 
                                    .Object@reader, PACKAGE="rbamtools")
        return(.Object)
})


bamReader <- function(filename, indexname, idx=FALSE, verbose=0)
{
    if(!is.logical(idx))
        stop("[bamReader] idx must be logical!")
    
    if(!is.numeric(verbose))
        stop("[bamReader] verbose must be numeric!")
    
    # Basic assumption is that only one file can be opened at once
    filename <- filename[1]
    idx <- idx[1]
    
    reader <- new("bamReader", filename)
    
    if((!idx) && missing(indexname))
    {
        if(verbose[1]==1)
            cat("[bamReader] Opened file '", basename(filename), "'.\n", sep="")
        else if(verbose[1]==2)
            cat("[bamReader] Opened file '", filename, "'.\n", sep="")
        return(reader)
    }
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    #  use indexname or set default
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    
    if(missing(indexname))
        idxfile <- paste(filename, "bai", sep=".")
    else
        idxfile <- indexname
    
    loadIndex(reader, idxfile)
    
    if(verbose[1]==1)
    {
        cat("[bamReader] Opened file '", basename(filename),
                    "' and index '", basename(idxfile), "'.\n", sep="")
    }else if(verbose[1]==2)
    {
        cat("[bamReader] Opened file '", filename,
                    "' and index '", idxfile, "'.\n", sep="")
    }
    
    return(reader)
}

setMethod("filename", "bamReader", function(object) return(object@filename))

setMethod("isOpen", signature="bamReader", definition=function(con, rw="")
{
    return(!(.Call("is_nil_externalptr", con@reader, PACKAGE="rbamtools")))
})

setMethod(f="bamClose", signature="bamReader", definition=function(object)
{
    if(!.Call("is_nil_externalptr", object@index, PACKAGE="rbamtools"))
    {
        .Call("bam_reader_unload_index", object@index, PACKAGE="rbamtools")
    }
    invisible(.Call("bam_reader_close", object@reader, PACKAGE="rbamtools"))
})


setMethod("show","bamReader", function(object)
{
    bm <- Sys.localeconv()[7]
    w <- 20
    r <- "right"
    cat("Class       : ", format(class(object)  , w=w, j=r)                       , "\n", sep="")
    cat("Filename    : ", format(basename(object@filename), w=w, j=r)             , "\n", sep="")
    cat("File status : ", format(ifelse(isOpen(object), "Open", "Closed"), w=w, j=r), "\n", sep="")
    cat("Index status: ", format(ifelse(indexInitialized(object), "Initialized", "Not initialized"), w=w, j=r), "\n", sep="")
    
    if(isOpen(object))
    {
        cat("RefCount    : ", format(getRefCount(object), w=w, j=r)                  , "\n", sep="")
    }
    return(invisible())
})

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#   End: Opening and closing a BAM-File for reading
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#
#   Header related functions
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#   This is one standard Method for creation of bamHeader
#   and is used as a simple way to pass a header to a new
#   instance of bamWriter
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

setMethod(f="getHeader", signature="bamReader", definition=function(object)
{
    if(!isOpen(object))
        stop("reader must be opened! Check with 'isOpen(reader)'!")
    
    return(new("bamHeader", .Call("bam_reader_get_header", object@reader))) 
})


setMethod(f="getHeaderText", signature="bamReader", definition=function(object)
{
    if(!isOpen(object))
        stop("reader must be opened! Check 'isOpen(reader)'!")
    
    return(new("bamHeaderText", .Call("bam_reader_get_header_text",
                                    object@reader, PACKAGE="rbamtools")))
})

#  + + + + + + + + + + + + + + + + + # 
#  getRefCount
#  + + + + + + + + + + + + + + + + + # 
setMethod(f="getRefCount", signature="bamReader",
                                            definition=function(object)
{
    if(!isOpen(object))
        stop("reader must be opened! Check with 'isOpen(reader)'!")
    
    return(.Call("bam_reader_get_ref_count",
                                object@reader, PACKAGE="rbamtools"))
})


#  + + + + + + + + + + + + + + + + + # 
#  getRefData
#  + + + + + + + + + + + + + + + + + # 
setMethod(f="getRefData", signature="bamReader", definition=function(object)
{
    if(!isOpen(object))
        stop("reader must be opened! Check with 'isOpen(reader)'!")
    
    return(.Call("bam_reader_get_ref_data", 
                                object@reader, PACKAGE="rbamtools"))
})


#  + + + + + + + + + + + + + + + + + # 
#  getRefCoords: Returns coordinates 
#  of entire reference for usage with
#  bamRange, gapList or siteList
#  function.
#  + + + + + + + + + + + + + + + + + # 
setMethod(f="getRefCoords",
                        signature="bamReader", definition=function(object,sn)
{
    if(!is.character(sn))
        stop("sn must be character!")
    
    if(length(sn) > 1)
        stop("sn must have length 1!")
    
    ref <- getRefData(object)
    id <- which(sn==ref$SN)
    
    if(length(id) == 0)
        stop("No match for sn in ref-data-SN!")
    
    coords <- c(ref$ID[id], 0, ref$LN[id])
    names(coords) <- c("refid", "start", "stop")
    return(c(ref$ID[id], 0, ref$LN[id]))
})

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#   End Header related functions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#   Index related functions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

setMethod(f="createIndex", signature="bamReader",
    definition=function(object,idx_filename)
{
    if(missing(idx_filename))
        idx_filename <- paste(object@filename, ".bai", sep="")
    
    return(invisible(.Call("bam_reader_create_index",
                            path.expand(object@filename),
                            path.expand(idx_filename), PACKAGE="rbamtools")))
})

setMethod("loadIndex", signature="bamReader", 
                                        definition=function(object,filename)
{
    if(!is.character(filename))
        stop("Filename must be character!\n")
    
    if(!file.exists(filename))
        stop("Index file \"", filename, "\" does not exist!\n")
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    #  Set index Variable in given bamReader object:
    #  Read object name, create expression string and evaluate in parent frame
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    reader <- deparse(substitute(object))
    
    extxt <- paste(reader,"@index<-.Call(\"bam_reader_load_index\",\"",
                path.expand(filename), "\",PACKAGE=\"rbamtools\")", sep="")
    
    eval.parent(parse(text=extxt))
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    #  Return true if bamReader@index!=NULL (parent frame)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    extxt <- paste(".Call(\"is_nil_externalptr\",", reader,
                                    "@index,PACKAGE=\"rbamtools\")", sep="")
    
    return(invisible(!eval.parent(parse(text=extxt))))
})


setMethod("indexInitialized", signature="bamReader", definition=function(object)
{ return(!(.Call("is_nil_externalptr", object@index, PACKAGE="rbamtools"))) })



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#  Deprecated functions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

setMethod("create.index", "bamReader", function(object, idx_filename)
{
    message("[create.index] Will soon be deprecated. Use createIndex.")
    return(createIndex(object, idx_filename))
    #.Deprecated(new="createIndex",package="rbamtools")
})


setMethod("load.index", "bamReader", function(object, filename)
{
    message("[load.index] Will soon be deprecated. Use loadIndex.")
    #.Deprecated(new="loadIndex", package="rbamtools")
    
    if(!is.character(filename))
        stop("Filename must be character!\n")
    
    if(!file.exists(filename))
        stop("Index file \"", filename, "\" does not exist!\n")
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    #  Set index Variable in given bamReader object:
    #  Read object name, create expression string and evaluate in parent frame
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    reader <- deparse(substitute(object))
    
    extxt <- paste(reader,"@index<-.Call(\"bam_reader_load_index\",\"",
                   path.expand(filename), "\",PACKAGE=\"rbamtools\")", sep="")
    
    eval.parent(parse(text=extxt))
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    #  Return true if bamReader@index!=NULL (parent frame)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    extxt <- paste(".Call(\"is_nil_externalptr\",", reader,
                   "@index,PACKAGE=\"rbamtools\")", sep="")
    
    return(invisible(!eval.parent(parse(text=extxt))))
})


setMethod("index.initialized", "bamReader", function(object)
{
    message("[index.initialized] Will soon be deprecated. Use indexInitialized")
    #.Deprecated(new="indexInitialized", package="rbamtools")
    return(indexInitialized(object))
})

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

setMethod(f="bamSort", signature="bamReader",
        definition=function(object, prefix="sorted", byName=FALSE,
                    maxmem=1e+9, path=dirname(filename(object)))
{
    if(!isOpen(object))
        stop("bamReader must be opened!")
    
    if(!is.logical(byName))
        stop("[bamSort] byName must be logical!")
    
    if(length(byName) > 1)
        stop("[bamSort] byName must have length 1!")
    
    if(!is.numeric(maxmem))
        stop("[bamSort] maxmem must be numeric!")
    
    if(length(maxmem)>1)
        stop("[bamSort] maxmem must have length 1!")
            
    maxmem <- floor(maxmem)
    message("[bamSort] Filename: ", object@filename)
    message("[bamSort] Prefix  : ", prefix)
    message("[bamSort] Maxmem  : ", maxmem)
    message("[bamSort] By Name : ", byName)
    
    .Call("bam_reader_sort_file",
            object@filename,
            path.expand(file.path(path,prefix)),
            maxmem, byName, PACKAGE="rbamtools")
    
    cat("[bamSort] Sorting finished.\n")
    
    return(invisible(paste(prefix, ".bam", sep="")))
})

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#   End Index related functions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #


# getNextAlign
setMethod(f="getNextAlign", signature="bamReader", definition=function(object)
{
    ans <- .Call("bam_reader_get_next_align", 
                                        object@reader, PACKAGE="rbamtools")
    
    if(is.null(ans))
        return(invisible(NULL))
    else
        return(new("bamAlign", ans))
})



setMethod("reader2fastq", "bamReader", 
        function(object, filename, which, append=FALSE)
{
    message("[reader2fastq] Will soon be deprecated. Use readerToFastq")
    readerToFastq(object, filename, which, append)
})

setMethod("readerToFastq", "bamReader", 
          function(object, filename, which, append=FALSE)
{
    if(!isOpen(object))
        stop("Reader must be opened!")
    
    if(!is.logical(append))
        stop("'append' must be logical!")
    
    if(!is.character(filename))
        stop("'filename' must be character!")
              
    if(missing(which))
    {
        return(invisible(.Call("bam_reader_write_fastq", object@reader,
            filename, append, PACKAGE="rbamtools")))
    }else{
        if(!is.numeric(which))
            stop("'which' argument must be numeric!")
        
        ans <- .Call("bam_reader_write_fastq_index", object@reader, filename,
            as.integer(sort(unique(which))), append, PACKAGE="rbamtools")
        
        if(ans < length(which))
            message("[readerToFastq] EOF reached.\n")
        
        message("[readerToFastq]", ans, "records written.\n")
            return(invisible(ans))
    }
})


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#  Reading gap-lists
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
setMethod("gapList", "bamReader", function(object, coords)
{
    if(!indexInitialized(object))
        stop("Reader must have initialized index!")
    
    return(new("gapList", object, coords))
})

setMethod("siteList", "bamReader", function(object, coords)
{
    if(!indexInitialized(object))
        stop("Reader must have initialized index!")
    
    return(new("gapSiteList", object, coords))
})

setMethod("bamGapList", "bamReader", function(object)
{
    if(!indexInitialized(object))
        stop("Reader must have initialized index!")
    
    return(new("bamGapList", object))
})


setMethod("rewind", "bamReader", function(object)
{
    return(invisible(.Call("bam_reader_seek", 
                    object@reader, object@startpos, PACKAGE="rbamtools")))
})


setMethod("bamSave", "bamReader", function(object, writer)
{
    if(!is(writer, "bamWriter"))
        stop("'writer' must be 'bamWriter'!")
    
    if(!isOpen(object))
        stop("'reader' is not open! Check 'isOpen'!")
    
    if(!isOpen(writer))
        stop("'writer' is not open! Check 'isOpen'!")
    
    # Saving old reading position
    oldpos <- .Call("bam_reader_tell", object@reader, PACKAGE="rbamtools")
    
    # Reset reader to start position
    .Call("bam_reader_seek",
                    object@reader, object@startpos, PACKAGE="rbamtools")
    
    nAligns <- .Call("bam_reader_save_aligns",
                    object@reader, writer@writer, PACKAGE="rbamtools")
    
    bm <- Sys.localeconv()[7]
    
    message("[bamSave.bamReader] Saving ", format(nAligns, big.mark=bm),
            " to file '", basename(writer@filename), "' finished.\n", sep="")
    
    .Call("bam_reader_seek", object@reader, oldpos, PACKAGE="rbamtools")
    
    return(invisible(nAligns))
})


setMethod("bamCopy", "bamReader", function(object, writer, refids, verbose=FALSE)
{
    if(!is(writer, "bamWriter"))
        stop("writer must be 'bamWriter'!")
    
    if(!isOpen(object))
        stop("reader is not open! Check 'isOpen'!")
    
    if(!isOpen(writer))
        stop("writer is not open! Check 'isOpen'!")
    
    if(!indexInitialized(object))
        stop("reader must have initialized index! Check 'indexInitialized'!")
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    #  Check refids argument: When missing copy all ref's
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    ref <- getRefData(object)
    if(missing(refids))
    {
        refids <- ref$ID
        n <- length(refids)
        mtc <- 1:n
    }
    else
    {
        mtc <- match(refids, ref$ID)
        if(any(is.na(mtc)))
            stop("refids must be subset of Reference-ID's! Check 'getRefData'!")
        n <- length(refids)    
    }
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    #  Copy aligns with bamRanges as intermediate buffer
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    bm <- Sys.localeconv()[7]
    nAligns <- 0
    for(i in 1:n)
    {
        range <- bamRange(object, 
                        c(ref$ID[mtc[i]], 0, ref$LN[mtc[i]]), complex=FALSE)
        
        nAligns <- nAligns+size(range)
        if(verbose)
        {
            message("[bamCopy.bamReader] i: ", i, "\tCopying ", 
            format(size(range), big.mark=bm, width=10), 
                    " aligns for Reference '", ref$SN[mtc[i]], "'.\n", sep="")
        }
        
        bamSave(writer, range, ref$ID[mtc[i]])
        rm(range)
        gc()
    }
    message("[bamCopy.bamReader] Copying ", 
        format(nAligns, big.mark=bm, width=10), " aligns finished.\n", sep="")
})


setMethod("extractRanges", "bamReader",
    definition=function(object, ranges, filename, complex=FALSE, header, idxname)
{
    if(!isOpen(object))  
        stop("Provided reader must be opened!")
    
    if(!indexInitialized(object))
        stop("Provided reader must have initialized index!")
    
    if(missing(header))
    {
        header <- getHeader(object)
    }else{
        if(!is(header, "bamHeader"))
            stop("[extractRanges] header must be of class 'bamHeader'")
        
        message("[extractRanges] bamHeader provided. Slot 'headerLine' will be changed (SO: unknown). Slot 'refSeqDict' will be overwritten.\n")
    }
    if(!is.data.frame(ranges))
        stop("[extractRanges] ranges must be 'data.frame'!")
    
    if(!is.logical(complex))
        stop("[extractRanges] complex must be logical!")
    
    if(length(complex) > 1)
        stop("[extractRanges] complex must have length 1!")
    
    # Preparing ranges table
    if(!all(is.element(c("seqid", "start", "end"), names(ranges))))
        stop("ranges argument must contain columns 'seqid', 'start', 'end'!")
    
    
    # Preparing filenames
    file_prefix <- sub("^([^.]*).*", "\\1", basename(filename))
    unsort_filename <- file.path(dirname(filename),
                paste("unsort", paste(file_prefix, "bam", sep="."), sep="_"))
    
    filename <- file.path(dirname(filename), paste(file_prefix, "bam", sep="."))
    
    message("[extractRanges] Provided filename is changed to '", filename,
                                        "' (see help for 'bamSort').", sep="")
    
    if(missing(idxname))
    {
        idxname <- paste(filename, "bai", sep=".")
    }else{
        if(!is.character(idxname))
            stop("[extractRanges] idxname must be character!\n")
    }
    bm <- Sys.localeconv()[7]
    

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    #  prepare range data
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    # It is essential for indexing, that all ref-ID's which
    # occur in aligns are also (implicitly) present in the
    # reference sequence dictionary (RSD) section.
    #
    # E.g. when there is an align which has refid 4, there
    # must be at least 5 entries in RSD because they are
    # indexed implicitly (that is: there is no entry in RSD
    # which says refid=4).
    #
    # Their ID is identified with the numbers
    # 0 to [(number of Entries in RSD)-1].
    #
    # Otherwise samtools indexing crashes without warning.
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    
    gp <- ranges[, c("seqid", "start", "end")] 
    rd <- getRefData(object)
    mtc <- match(gp$seqid, rd$SN)
    isna <- is.na(mtc)
    if(all(isna))
    {
        message("[extractRanges] No matching seqids for genes:")
        print(gp)
        message("[extractRanges] No output generated.\n")
        return(invisible())
    }
    
    if(any(isna))
    {
        message("[extractRanges] Missing seqid matches. Skipping following genes:")
        print(gp[isna, ])
        gp <- gp[!isna, ]
    }
    
    n <- dim(gp)[1]
    gp$old_ID <- rd$ID[mtc]
    gp$LN <- rd$LN[mtc]
    
    # Provide (unique) new ID's
    renew <- data.frame(old=sort(unique(gp$old_ID)))
    nid <- dim(renew)[1]
    renew$new <- 0:(nid-1)
    mtc <- match(gp$old_ID, renew$old)
    
    gp$new_ID <- renew$new[mtc]
    
    new_rd <- merge(renew, rd, by.x="old", by.y="ID")
    nref <- dim(new_rd)[1]
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    #  create new header
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    
    htxt <- getHeaderText(header)
    hl <- headerLine(htxt)
    setVal(hl, "SO", "unsorted")
    
    rsd <- new("refSeqDict")
    for(i in 1:nref)
        addSeq(rsd, SN=new_rd$SN[i], LN=new_rd$LN[i])
    headerLine(htxt) <- hl
    refSeqDict(htxt) <- rsd
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    #  Writing aligns
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    
    cat("[extractRanges] Writing aligns to temporary file '", 
                                    unsort_filename, "'.\n", sep="")
    
    writer <- bamWriter(bamHeader(htxt), unsort_filename)
    
    nAligns <- 0
    for(i in 1:n)
    {
        range <- bamRange(object, c(gp$old_ID[i], gp$start[i], gp$end[i]))
        if(size(range)==0)
        {
            message("No aligns found for gene '", gp$gene_name[i], "'.")
        }else{
            bamSave(writer, range, refid=gp$new_ID[i])
            nAligns <- nAligns+size(range)
        }
    }
    message("[extractRanges] Writing of", format(nAligns, big.mark=bm),
                                                        "aligns finished.")
    
    bamClose(writer)
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    #  Sorting BAM output file
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    
    cat("[extractRanges] Sorting:\n")
    nread <- bamReader(unsort_filename)
    if(!isOpen(nread))
        stop("unsorted bam file not found!")
    bamSort(nread, prefix=file_prefix)
    bamClose(nread)
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    #  Creating index for ouput file
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    message("[extractRanges] Creating index '", basename(idxname), "'.", sep="")
    
    nread <- bamReader(filename)
    createIndex(nread, idx_filename=idxname)
    message("[extractRanges] Finished.\n")
    message("[extractRanges] You may want to delete file '", 
                                    basename(unsort_filename), "'.\n", sep="")
    return(invisible(nAligns))
})



setMethod("bamCount", signature="bamReader", definition=function(object, coords)
{
    if(!indexInitialized(object))
        stop("[bamCount] reader must have initialized index! Use 'loadIndex'!")
    
    if(missing(coords))
        stop("[bamCount] coords is not optional!")
    
    if(!is.numeric(coords))
        stop("[bamCount] coords must be numeric")
    
    res <- .Call("bam_count", object@reader, 
                                    object@index, coords, PACKAGE="rbamtools")
    
    names(res) <- c("M", "I", "D", "N", "S", "H", "P", "=", "X", "nAligns")
    return(res)
})


setMethod("bamCountAll", "bamReader", function(object, verbose=FALSE)
{
    if(!isOpen(object))
        stop("reader is not open! Check 'isOpen'!")
    
    if(!indexInitialized(object))
        stop("reader must have initialized index! Check 'indexInitialized'!")
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    #  Check refids argument: When missing copy all ref's
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    ref <- getRefData(object)
    nr <- nrow(ref)
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    #  Count first refid
    #  and read size and names of result
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    if(verbose)
        cat("[bamCountAll] Counting ", ref$SN[1],  "\t[ 1/", nr, "]", sep="")
    count <- bamCount(object, c(ref$ID[1], 0, ref$LN[1]))
    nc <- length(count)
    
    mtx <- matrix(numeric(nc*nr), ncol=nc)
    colnames(mtx) <- names(count)
    rownames(mtx) <- ref$SN
    mtx[1, ] <- count
    
    if(nr > 1)
    {
        for(i in 2:nr)
        {
            if(verbose)
            {
                cat("\r[bamCountAll] Counting ", ref$SN[i], 
                            "\t[", format(i, width=2), "/", nr, "]", sep="")
            }
            mtx[i, ] <- bamCount(object, c(ref$ID[i], 0, ref$LN[i]))
        }
    }
    
    if(verbose)
        cat("\n[bamCountAll] Finished.\n")
    res <- as.data.frame(mtx)
    res$ID <- ref$ID
    res$LN <- ref$LN
    return(res)
})


setMethod("nucStats", "bamRange", function(object)
{
    m <- countNucs(object)
    gcc <- gc_content(m[1], m[2], m[3], m[4])
    at_gc <- at_gc_ratio(m[1], m[2], m[3], m[4])
    dfr<-data.frame(
                        nAligns=size(object),
                        A=m[1],
                        C=m[2],
                        G=m[3],
                        T=m[4],
                        N=m[5],
                        gcc=gcc,
                        at_gc_ratio=at_gc
                    )
    
    refname <- .Call("bam_range_get_refname",
                        object@range, PACKAGE="rbamtools")
    
    if(!is.null(refname))
        row.names(dfr) <- refname
    
    return(dfr)
})



setMethod("nucStats", "bamReader", function(object)
{
    if(!isOpen(object))
        stop("Reader must be open (check 'isOpen')!")
    if(!indexInitialized(object))
        stop("Reader must have initialized index (use 'loadIndex')!")
    
    
    ref <- getRefData(object)
    n <- nrow(ref)
    m <- matrix(0, nrow=n, ncol=5)
    nAligns <- numeric(n)
    for(i in 1:n)
    {
        range <- bamRange(object, c(ref$ID[i], 0, ref$LN[i]))
        nAligns[i] <- size(range)
        m[i, ] <- countNucs(range)
    }
    dfr <- data.frame(nAligns=nAligns, 
                            A=m[, 1], C=m[, 2], G=m[, 3], T=m[, 4], N=m[, 5])
    
    dfr$gcc <- gc_content(dfr$A, dfr$C, dfr$G, dfr$T)
    
    dfr$at_gc_ratio <- at_gc_ratio(dfr$A, dfr$C, dfr$G, dfr$T)
    
    rownames(dfr) <- ref$SN
    return(dfr)
})

setMethod("nucStats", "character", 
        definition=function(object, idxInfiles=paste(object, ".bai", sep=""))
{
    if(any(!file.exists(object)))
        stop("[nucStats] Files (object) not found!")
    
    if(!is.character(idxInfiles))
        stop("[nucStats] idxInfiles must be character!")
    
    if(any(!file.exists(idxInfiles)))
        stop("[nucStats] Files (idxInfiles) not found!")
    
    n <- length(object)
    if(length(idxInfiles)!=n)
        stop("[nucStats] infiles and idxInfiles must have same length!")
    
    res <- data.frame(nAligns=numeric(n), A=numeric(n),
                C=numeric(n), G=numeric(n), T=numeric(n), N=numeric(n))
    
    for(i in 1:n)
    {
        reader <- bamReader(object[i])
        loadIndex(reader, idxInfiles[i])
        nc <- nucStats(reader)
        res[i, ] <- lapply(nc[, 1:6], sum)
    }
    res$gcc <- gc_content(res$A, res$C, res$G, res$T)
    res$at_gc_ratio <- at_gc_ratio(res$A, res$C, res$G, res$T)
    rownames(res) <- 1:n
    return(res)
})

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#
#  bamHeader
#  Description: See SAM File Format Specification (v1.4-r985)
#  September 7, 2011, Section 1.3
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #


setMethod("initialize", "bamHeader", function(.Object, extptr)
{
    if(!is(extptr, "externalptr"))
        stop("extptr must be externalptr!")
    .Object@header <- extptr
    return(.Object)
})

setMethod(f="getHeaderText", signature="bamHeader", definition=function(object)
{
    return(new("bamHeaderText", .Call("bam_header_get_header_text", 
                                object@header, PACKAGE="rbamtools"))) })

setMethod("as.character", "bamHeader", function(x, ...)
{
    .Call("bam_header_get_header_text", x@header, PACKAGE="rbamtools")
})

setMethod("show", "bamHeader", function(object)
{
    cat("An object of class \"", class(object), "\"\n", sep="")
    ht <- getHeaderText(object)
    hl <- headerLine(ht)
    dc <- refSeqDict(ht)
    
    cat("headerLine:\n")
    cat("VN:", hl@VN, "\n")
    cat("SO:", hl@SO, "\n")
    
    nsq <- length(dc@SN)
    cat("refSeqDict: (size ", nsq, ")\n", sep="")
    if(nsq>0)
    {
        cat("Seqs: ")
        for(i in 1:(pmin(nsq, 3)))
            cat(dc@SN[i], ",  ")
        cat("...\n")    
    }
})

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#  This is the main function for creating an instance of bamWriter
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

setMethod("bamWriter", "bamHeader", function(x, filename)
{
    if(!is.character(filename))
        stop("[bamWriter.bamHeader] filename must be character!")
    return(new("bamWriter", x, filename))
})


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# headerLine: Represents two entries: 
# Format version (VN) and sorting order(SO)
# Valid format for VN : /^[0-9]+\.[0-9]+$/.
# Valid entries for SO: unknown (default),  unsorted,  queryname,  coordinate.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

setMethod(f="initialize", signature="headerLine", 
                            definition=function(.Object, hl="", delim="\t")
{
    # Parses header line from header section
    if(!is.character(hl))
        stop("[headerLine.initialize] Argument must be string.\n")
    
    # Default object content (hl="" or character(0))
    if((length(hl)==1 && nchar(hl)==0) || length(hl)==0)
    {
        .Object@VN="1.4"
        .Object@SO="unknown"
        return(.Object)
    }
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    #   Split input string into tags
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    tags <- unlist(strsplit(hl, delim))
    
    #  Three tags!
    if(length(tags)!=3)
        stop("hl must contain three tags separated by '", delim, "'!\n")
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    #   First  tag: '@HD'
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    if(tags[1]!="@HD")
        stop("First tag of string must be @HD!\n")
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    #   Second tag: 'VN'
    #   TODO: Check Accepted format: /^[0-9]+\.[0-9]+$/.  
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    if(substr(tags[2], 1, 2)!="VN")
        stop("Second tag of string must be VN!\n")
    
    .Object@VN=substring(tags[2], 4)
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    #  Third   tag: 'SO'
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    if(substr(tags[3], 1, 2)!="SO")
        stop("Third tag of string must be SO!\n")
    
    str <- substring(tags[3], 4)
    if(str=="coordinate")
        .Object@SO <- "coordinate"
    else if(str=="unknown")
        .Object@SO <- "unknown"
    else if(str=="unsorted")
        .Object@SO <- "unsorted"
    else if(str=="queryname")
        .Object@SO <- "queryname"
    
    return(.Object)
})

setMethod("getHeaderText", "headerLine", function(object, delim="\t")
    {return(paste("@HD\tVN:", object@VN, "\tSO:", object@SO, sep=""))})

setMethod("getVal", signature="headerLine", definition=function(object, member)
{
    if(!is.character(member))
        stop("[getVal.headerLine] Member must be character!\n")
    if(member=="VN")
        return(object@VN)
    if(member=="SO")
        return(object@SO)
    stop("Member '", member, "' must be 'VN' or 'SO'!\n")
})

setMethod("setVal", signature="headerLine",
                                definition=function(object, members, values)
{
    if(!is.character(members) || !is.character(values))
        stop("Members and values must be character!\n")
    if(length(members)!=length(values))
        stop("Members and values must have same length!\n")
    
    tagLabs <- c("VN", "SO")
    mtc <- match(members, tagLabs)
    if(any(is.na(mtc)))
        stop("Member names must be valid Header line entries!\n")
    
    n <- length(members)
    if(n>2)
        stop("Only two members can be set!\n")
    
    obj <- deparse(substitute(object))
    for(i in 1:n)
    {
        txt <- paste(obj, "@", members[i], " <- '", values[i], "'", sep="")
        eval.parent(parse(text=txt))
    }
    return(invisible())
})

setMethod("as.list", signature="headerLine", definition=function(x, ...)
    {return(list(VN=x@VN, SO=x@SO))})

setMethod("show", "headerLine", function(object)
{
    cat("An object of class \"", class(object), "\"\n", sep="")
    cat("VN: ", object@VN, "\nSO: ", object@SO, "\n", sep="")
})

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#  End headerLine
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#   refSeqDict: Reference Sequence Dictionary
#   Represents a variable number of Ref Seqs
#   Valid Members (Entries for each sequence, stored in a data.frame):
#   SN Reference sequence name
#   LN Reference sequence length
#   AS Genome assembly identifier
#   M5 MD5 checksum of the sequence
#   SP Species
#   UR URI of the sequence
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

setMethod(f="initialize", signature="refSeqDict", 
                            definition=function(.Object, hsq="", delim="\t")
{
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    #  Parses Reference sequence dictionary of header-text
    #  hsq= Vector of characters,  each representing one Ref-Sequence
    #  length(hsq) = number of Ref-Sequences
    #  Each Ref-string contains 'internally' [tab] delimited seqments:
    #                "SN:ab\tLN:12\tAS:ab\tM5:12\tSP:ab\tUR:ab"
    #  It's allowed to skip segments
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    
    if(!is.character(hsq))
        stop("[refSeqDict.initialize] hsq must be character!")
    
    n <- length(hsq)
    # Return empty object when no input string is given
    if((n==1 && nchar(hsq)==0) || n==0)
        return(.Object)
    
    .Object@SN <- character(n)
    .Object@LN <- numeric(n)
    .Object@AS <- character(n)
    .Object@M5 <- numeric(n)
    .Object@SP <- character(n)
    .Object@UR <- character(n)
    
    labels <- c("SN", "LN", "AS", "M5", "SP", "UR")
    for(i in 1:n)
    {
        # Containes separated tags for one sequence
        seq <- unlist(strsplit(hsq[i], delim))
        if(seq[1]!="@SQ")
        stop("First segment in Ref-sequence tag must be '@SQ'!")
        seq <- seq[-1]
        
        # Contains column number in dict@df for each tag
        cols <- match(substr(seq, 1, 2), labels)
        m <- length(cols)
        
        for(j in 1:m)
        {
            txt <- substr(seq[j], 4, nchar(seq[j]))
            # Empty entries are skipped (to avoid errors)
            if(nchar(txt)>0)
            {
                if(cols[j]==1)
                    .Object@SN[i] <- txt
                else if(cols[j]==2)
                {
                    # Try to convert into numeric value
                    numb <- suppressWarnings(as.numeric(txt))  
                    if(is.na(numb))
                    {
                        warning("[refSeqDict.initialize] No numeric value for LN: '",
                                                txt, "'!\n", sep="")
                    }else{
                        .Object@LN[i] <- numb
                    }
                }
                else if(cols[j]==3)
                    .Object@AS[i] <- txt
                else if(cols[j]==4)
                {
                    # Try to convert into numeric value
                    numb <- suppressWarnings(as.numeric(txt))  
                    if(is.na(numb))
                    {
                        warning("[refSeqDict.initialize] No numeric value for LN: '",
                                                txt, "'!\n", sep="")
                    }else{
                        .Object@M5 <- numb
                    }
                }
                else if(cols[j]==5)
                    .Object@SP[i] <- txt
                else if(cols[j]==6)
                    .Object@UR[i] <- txt
            }
        }
    }
    return(.Object)
})

setMethod(f= "[", signature="refSeqDict", definition=function(x, i)
{
    rsd <- new("refSeqDict")
    rsd@SN <- x@SN[i]
    rsd@LN <- x@LN[i]
    rsd@AS <- x@AS[i]
    rsd@M5 <- x@M5[i]
    rsd@SP <- x@SP[i]
    rsd@UR <- x@UR[i]
    return(rsd)
})

setMethod(f="dim", signature="refSeqDict", definition=function(x)
                            {return(c(length(x@SN), 6))})


setMethod("removeSeqs", signature="refSeqDict", definition=function(x, rows)
{
    #  Removes given rows (=Sequences) from Dictionary
    #  so they are excluded from header
    n <- length(x@SN)
    if(!is.numeric(rows))  
        stop("[removeSeqs.refSeqDict] Sequence indices must be numeric!")
    rows <- as.integer(rows)
    
    if(any(rows)<1)
        stop("[removeSeqs.refSeqDict] Sequence indices must be positive!")
    if(any(rows)>n)
        stop("[removeSeqs.refSeqDict] Sequence indices must be <", n, "!")
    
    # Execute per eval in parent.frame
    if(length(rows)>1)
        rmv <- paste("c(", paste(rows, collapse=", "), ")", sep="")
    else
    rmv <- rows
    
    obj <- deparse(substitute(x))
    dictcol <- paste(obj, "@SN", sep="")
    eval.parent(parse(text=paste(dictcol, "<-", dictcol, "[-", rmv, "]", sep="")))
    
    dictcol<-paste(obj, "@LN", sep="")
    eval.parent(parse(text=paste(dictcol, "<-", dictcol, "[-", rmv, "]", sep="")))
    
    dictcol<-paste(obj, "@AS", sep="")
    eval.parent(parse(text=paste(dictcol, "<-", dictcol, "[-", rmv, "]", sep="")))
    
    dictcol<-paste(obj, "@M5", sep="")
    eval.parent(parse(text=paste(dictcol, "<-", dictcol, "[-", rmv, "]", sep="")))
    
    dictcol<-paste(obj, "@SP", sep="")
    eval.parent(parse(text=paste(dictcol, "<-", dictcol, "[-", rmv, "]", sep="")))
    
    dictcol<-paste(obj, "@UR", sep="")
    eval.parent(parse(text=paste(dictcol, "<-", dictcol, "[-", rmv, "]", sep="")))
    return(invisible())
})


setMethod("addSeq", signature="refSeqDict", 
                definition=function(object, SN, LN, AS="", M5=0, SP="", UR="")
{
    index <- length(object@SN)+1
    obj <- deparse(substitute(object))
    colidx <- paste("[", index, "]", sep="")
    
    # Appends new Sequence (row) at the end
    dictcol <- paste(obj, "@SN", colidx, sep="")
    eval.parent(parse(text=paste(dictcol, "<-'", SN, "'", sep="")))
    
    dictcol <- paste(obj, "@LN", colidx, sep="") 
    eval.parent(parse(text=paste(dictcol, "<-", LN, sep="")))
    
    dictcol <- paste(obj, "@AS", colidx, sep="") 
    eval.parent(parse(text=paste(dictcol, "<-'", AS, "'", sep="")))
    
    dictcol <- paste(obj, "@M5", colidx, sep="") 
    eval.parent(parse(text=paste(dictcol, "<-", M5, sep="")))
    
    dictcol <- paste(obj, "@SP", colidx, sep="") 
    eval.parent(parse(text=paste(dictcol, "<-'", SP, "'", sep="")))
    
    dictcol <- paste(obj, "@UR", colidx, sep="") 
    eval.parent(parse(text=paste(dictcol, "<-'", UR, "'", sep="")))
    
    return(invisible())
})

setMethod("getHeaderText", signature="refSeqDict", 
                                    definition=function(object, delim="\t")
{
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    #  Returns Ref Data String (can be used for creating new BAM 
    #  file via bamWriter)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
                                        
    labels <- c("SN", "LN", "AS", "M5", "SP", "UR")
    n <- length(object@SN)
    
    if(n==0)
        return(character(0))
    
    seqs <- character(n)
    
    for(i in 1:n)
    {
        ans <- "@SQ"    
        if(nchar(object@SN[i])>0)
            ans <- paste(ans, delim, "SN:", object@SN[i], sep="")
        if(object@LN[i]>0)
            ans <- paste(ans, delim, "LN:", object@LN[i], sep="")
        if(nchar(object@AS[i])>0)
            ans <- paste(ans, delim, "AS:", object@AS[i], sep="")
        if(object@M5[i]>0)
            ans <- paste(ans, delim, "M5:", object@M5[i], sep="")
        if(nchar(object@SP[i])>0)
            ans <- paste(ans, delim, "SP:", object@SP[i], sep="")
        if(nchar(object@UR[i])>0)
            ans <- paste(ans, delim, "UR:", object@UR[i], sep="")
        seqs[i] <- ans
    }
    return(paste(seqs, collapse="\n"))
})


#  Return first or last part of refSeqDict data.frame
#  S3 Generic is supplied via importFrom in NAMESPACE

setMethod("head", "refSeqDict", function(x, n=6L, ...)
{
    stopifnot(length(n) == 1L)
    if (n < 0L)
        stop("[head.refSeqDict] n<0!")
    
    m <- length(x@SN)
    if(m==0)
        cat("[head.refSeqDict] Empty object.\n")
    
    n <- min(n, m)
    if(n == 0L)
        return(as.data.frame(new("refSeqDict")))
    else
        return(as.data.frame(x)[1:n, ])
})

# S3 Generic is supplied via importFrom in NAMESPACE
setMethod("tail", "refSeqDict", definition=function(x, n=6L, ...)
{
    stopifnot(length(n) == 1L)
    if (n < 0L)
        stop("[tail.refSeqDict] n<0!")
    m <- length(x@SN)
    if(m==0)
        cat("[tail.refSeqDict] Empty object.\n") 
    n <- min(n, m)
    if(n == 0L)
        return(as.data.frame(new("refSeqDict")))
    else
    {
        n <- m-n+1
        return(x@df[n:m, ])
    }
})

setMethod("show", "refSeqDict", function(object)
{
    if(length(object@SN)>0)
    {
        cat("An object of class \"", class(object), "\"\n", sep="")
        print(head(object))
    }else{
        cat("An empty object of class \"", class(object), "\".\n", sep="")
    }
})

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#   End refSeqDict
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#  headerReadGroup
#  ReadGroup
#  ID Read Group identifier
#  CN Name of sequencing center
#  DS Description
#  FO Flow order
#  KS Nucleotides corresponding to key sequence of each read
#  LB Library
#  PG Programs used for processing the Read Group
#  PI Predicted median insert size
#  PL Sequencing Platform:
#     CAPILLARY, LS454, ILLUMINA, SOLID, HELICOS, IONTORRENT or PACBIO
#  SM Sample name.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #



setMethod(f="initialize", signature="headerReadGroup",  
                                definition=function(.Object, hrg="", delim="\t")
{
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    #  Parses Read-Group part of Header data. See Sam Format 
    #  Specificatioin 1.3 (Header Section)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    
    if(!is.character(hrg))
        stop("[headerReadGroup.initialize] Argument must be string.\n")
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    #  Samtools file format says: Unordered multiple @RG lines are allowed
    #  Each @RG segment comes as one string in hrg
    #  Number of @RG segments = length(hrg)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    #  Split hgr into multiple @RG fragments
    #  In effect, hrg can be either given as vector with length > 1
    #  or as single vector with @RG entries separated by '\n'
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    hrg <- unlist(strsplit(hrg, "\n"))
    
    .Object@nrg <- length(hrg)
    
    .Object@ntl <- 12L                      # number of tags
    .Object@ID <- character(.Object@nrg)
    .Object@CN <- character(.Object@nrg)
    .Object@DS <- character(.Object@nrg)
    .Object@DT <- character(.Object@nrg)
    .Object@FO <- character(.Object@nrg)
    .Object@KS <- character(.Object@nrg)
    .Object@LB <- character(.Object@nrg)
    .Object@PG <- character(.Object@nrg)
    .Object@PI <- character(.Object@nrg)
    .Object@PL <- character(.Object@nrg)
    .Object@PU <- character(.Object@nrg)
    .Object@SM <- character(.Object@nrg)
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    #  Allows for empty object
    #  hrg="" or character(0)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    if((length(hrg)==1 && nchar(hrg)==0) || length(hrg)==0)
    {
        .Object@nrg <- 0L
        return(.Object)
    }
    
    tagLabs <- c("ID", "CN", "DS", "DT", "FO", "KS", "LB", "PG",
                                                    "PI", "PL", "PU", "SM")
    
    for(i in 1:(.Object@nrg))
    {
        #  + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + # 
        #  Split string into fields
        #  + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + # 
        tags <- unlist(strsplit(hrg[i], delim))
        
        if(tags[1] != "@RG")
            stop("First item of string must be @RG!\n")
        
        #  TODO: Routine does not check for:
        #  'Each @RG line must have a unique ID.' (SAM file format)
        
        tags <- tags[-1]
        ntags <- length(tags)
        for(j in 1:ntags)
        {
            if(substr(tags[j], 1, 2) == "ID")
                .Object@ID[i] <- substring(tags[j], 4)
            
            if(substr(tags[j], 1, 2) == "CN")
                .Object@CN[i] <- substring(tags[j], 4)
            
            if(substr(tags[j], 1, 2) == "DS")
                .Object@DS[i] <- substring(tags[j], 4)
            
            if(substr(tags[j], 1, 2) == "DT")
                .Object@DT[i] <- substring(tags[j], 4)
            
            if(substr(tags[j], 1, 2) == "FO")
                .Object@FO[i] <- substring(tags[j], 4)
            
            if(substr(tags[j], 1, 2) == "KS")
                .Object@KS[i] <- substring(tags[j], 4)
            
            if(substr(tags[j], 1, 2) == "LB")
                .Object@LB[i] <- substring(tags[j], 4)
            
            if(substr(tags[j], 1, 2) == "PG")
                .Object@PG[i] <- substring(tags[j], 4)
            
            if(substr(tags[j], 1, 2) == "PI")
                .Object@PI[i] <- substring(tags[j], 4)
            
            if(substr(tags[j], 1, 2) == "PL")
                .Object@PL[i] <- substring(tags[j], 4)
            
            if(substr(tags[j], 1, 2) == "PU")
                .Object@PU[i] <- substring(tags[j], 4)
            
            if(substr(tags[j], 1, 2) == "SM")
                .Object@SM[i] <- substring(tags[j], 4)
        }
    }
    return(.Object)
})


setMethod("show", "headerReadGroup",  function(object)
{
    n <- object@nrg
    if(n > 0)
    {
        cat("An object of class \"", class(object), "\"\n", sep="")
        for(i in 1:n)
        {
            if(nchar(object@ID[i]) > 0)
                cat("ID:", object@ID[i], "\n")
            
            if(nchar(object@CN[i]) > 0)
                cat("CN:", object@CN[i], "\n")
            
            if(nchar(object@DS[i]) > 0)
                cat("DS:", object@DS[i], "\n")
            
            if(nchar(object@DT[i]) > 0)
                cat("DT:", object@DT[i], "\n")
            
            if(nchar(object@FO[i]) > 0)
                cat("FO:", object@FO[i], "\n")
            
            if(nchar(object@KS[i]) > 0)
                cat("KS:", object@KS[i], "\n")
            
            if(nchar(object@LB[i]) > 0)
                cat("LB:", object@LB[i], "\n")
            
            if(nchar(object@PG[i]) > 0)
                cat("PG:", object@PG[i], "\n")
            
            if(nchar(object@PI[i]) > 0)
                cat("PI:", object@PI[i], "\n")
            
            if(nchar(object@PL[i]) > 0)
                cat("PL:", object@PL[i], "\n")
            
            if(nchar(object@PU[i]) > 0)
                cat("PU:", object@PU[i], "\n")
            
            if(nchar(object@SM[i]) > 0)
                cat("SM:", object@SM[i], "\n")
            
            if(n > i)
                cat("\n")
        }
    }else{
        cat("An empty object of class \"", class(object), "\"\n", sep="")
    }
    return(invisible())
})

setMethod("getHeaderText", signature="headerReadGroup", 
                                        definition=function(object, delim="\t")
{
    n <- object@nrg
    if(n == 0)
        return(character(0))
    
    rgtxt <- character(n)
    for(i in 1:n)
    {
        #  ID should always be present
        txt <- paste(delim,"ID:",object@ID[i], sep="")
        
        if(nchar(object@CN[i]) > 0)
            txt <- paste(txt, delim,  "CN:", object@CN[i], sep="")
        
        if(nchar(object@DS[i]) > 0)
            txt <- paste(txt, delim, "DS:", object@DS[i], sep="")
        
        if(nchar(object@DT[i]) > 0)
            txt <- paste(txt, delim, "DT:", object@DT[i], sep="")
        
        if(nchar(object@FO[i]) > 0)
            txt <- paste(txt, delim, "FO:", object@FO[i], sep="")
        
        if(nchar(object@KS[i]) > 0)
            txt <- paste(txt, delim, "KS:", object@KS[i], sep="")
        
        if(nchar(object@LB[i]) > 0)
            txt <- paste(txt, delim, "LB:", object@LB[i], sep="")
        
        if(nchar(object@PG[i]) > 0)
            txt <- paste(txt, delim, "PG:", object@PG[i], sep="")
        
        if(nchar(object@PI[i]) > 0)
            txt <- paste(txt, delim, "PI:", object@PI[i], sep="")
        
        if(nchar(object@PL[i]) > 0)
            txt <- paste(txt, delim, "PL:", object@PL[i], sep="")
        
        if(nchar(object@PU[i]) > 0)
            txt <- paste(txt, delim, "PU:", object@PU[i], sep="")
        
        if(nchar(object@SM[i]) > 0)
            txt <- paste(txt, delim, "SM:", object@SM[i], sep="")
        
        # Remove last "\t"
        rgtxt[i] <- paste("@RG",txt,"\n",sep="")
    }
    return(paste(rgtxt, collapse=""))
})

setMethod("getVal", signature="headerReadGroup", 
                                        definition=function(object, member)
{
    if(!is.character(member))
        stop("Member must be character!\n")
    
    
    tagLabs <- c("ID", "CN", "DS", "DT", "FO", "KS", "LB", 
                                                "PG", "PI", "PL", "PU", "SM")
    
    mtc <- match(member, tagLabs)
    
    if(any(is.na(mtc)))
        stop("Invalid member name!\n")
    
    l <- list()
    if(object@nrg == 0)
        return(l)
    
    for(i in 1:length(member))
    {
        if(mtc[i] == 1)
            l$ID <- object@ID
        else if(mtc[i] == 2)
            l$CN <- object@CN
        else if(mtc[i] == 3)
            l$DS <- object@DS
        else if(mtc[i] == 4)
            l$DT <- object@DT
        else if(mtc[i] == 5)
            l$FO <- object@FO
        else if(mtc[i] == 6)
            l$KS <- object@KS
        else if(mtc[i] == 7)
            l$LB <- object@LB
        else if(mtc[i] == 8)
            l$PG <- object@PG
        else if(mtc[i] == 9)
            l$PI <- object@PI
        else if(mtc[i] == 10)
            l$PL <- object@PL
        else if(mtc[i] == 11)
            l$PU <- object@PU
        else if(mtc[i] == 12)
            l$SM <- object@SM
    }
    return(l)
})

setMethod("setVal", signature="headerReadGroup", 
                                    definition=function(object, members, values)
{
    if(!is.character(members))
        stop("Member name must be character!\n")
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    #  Check for valid values list
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    
    if(!is.list(values))
        stop("Values must be given as list object")
    
    if(length(members)!=length(values))
        stop("Members and values must have same length!\n")
    
    if(any(lapply(values,length)!=object@nrg))
        stop("Length of values must equal number of read groups!")
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    #  Check for valid members entries
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    tagLabs <- c("ID", "CN", "DS", "DT", "FO", "KS", "LB", "PG", "PI", 
                                                            "PL", "PU", "SM")
    
    mtc <- match(members, tagLabs)
    
    if(any(is.na(mtc)))
    {
        stop("Members must be valid Read Group Entries ",
                                        "(See SAM Format Specification 1.3!\n")
    }
    
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    #  Create insertion code as string and parse in parent environment
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    n <- length(members)
    obj <- deparse(substitute(object))
    for(i in 1:n)
    {
        for(j in 1:object@nrg)
        {
            txt <- paste(obj, "@", members[i], "[", j, "]", "<-'",
                                                values[[i]][j], "'", sep="")
            eval.parent(parse(text=txt))            
        }

    }
    return(invisible())
})

setMethod("as.list", signature="headerReadGroup", 
                                            definition=function(x, ...)
{
    l <- list()
    l$ID <- x@ID
    l$CN <- x@CN
    l$DS <- x@DS
    l$DT <- x@DT
    l$FO <- x@FO
    l$KS <- x@KS
    l$LB <- x@LB
    l$PG <- x@PG
    l$PI <- x@PI
    l$PL <- x@PL
    l$PU <- x@PU
    l$SM <- x@SM
    return(l)
})


setMethod("addReadGroup", signature="headerReadGroup", 
                                            definition=function(object, l)
{
    if(!is.list(l))
        stop("'l' must be list")
    
    tagLabs <- c("ID", "CN", "DS", "DT", "FO", "KS", "LB", "PG", "PI", 
                            "PL", "PU", "SM")
    
    mtc <- match(names(l), tagLabs)
    if(any(is.na(mtc)))
        stop("All list names must be valid read group tags.")
    
    mtc <- match(tagLabs, names(l))
    if(is.na(mtc[1]))
       stop("There must be an ID given for new read group")
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    #  Increase number of read groups by 1
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    n <- object@nrg
    object@nrg <- n + 1L
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    #  Insert ID tag
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    if(n == 0)
    {
        object@ID <- l[[mtc[1]]]
    }else{
        object@ID <- c(object@ID, l[[mtc[1]]])        
    }
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    #  Eventually insert other tags
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    ins <- function(o, i)
    {
        if(!is.na(mtc[i]))
        {
            if(n == 0)
                o <- l[[mtc[i]]]
            else if(length(o) > 0)
                o <- c(o, l[[mtc[i]]])
            else
                o <- c(rep("", n), l[[mtc[i]]])
        }else{
            if(n == 0)
                o <- ""
            else
                o <- c(o,"")
        }
        
        return(o)
    }
    
    i <- 2
    object@CN <- ins(object@CN, i)
    i <- i + 1
    object@DS <- ins(object@DS, i)
    i <- i + 1
    object@DT <- ins(object@DT, i)
    i <- i + 1
    object@FO <- ins(object@FO, i)
    i <- i + 1
    object@KS <- ins(object@KS, i)
    i <- i + 1
    object@LB <- ins(object@LB, i)
    i <- i + 1
    object@PG <- ins(object@PG, i)
    i <- i + 1
    object@PI <- ins(object@PI, i)
    i <- i + 1
    object@PL <- ins(object@PL, i)
    i <- i + 1
    object@PU <- ins(object@PU, i)
    i <- i + 1
    object@SM <- ins(object@SM, i)

    return(object)
})



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#   End headerReadGroup
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#  headerProgram
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #


setMethod(f="initialize", signature="headerProgram", 
                        definition=function(.Object, hp="", delim="\t")
{
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    #  Parses Program part of Header data.
    #  See Sam Format Specificatioin 1.3 (Header Section)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    .Object@l <- list()
    
    if(!is.character(hp))
        stop("[headerProgram.initialize] Argument must be string.\n")
    
    # hp="" or character(0)
    if((length(hp)==1 && nchar(hp)==0)||length(hp)==0)
        return(.Object)
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    #  Split string into fields
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    tags <- unlist(strsplit(hp, delim))
    if(tags[1]!="@PG")
        stop("[headerProgram.initialize] First item of string must be @PG!\n")
    
    tags <- tags[-1]
    tagLabs <- c("ID", "PN", "CL", "PP", "DS", "VN")
    n <- length(tags)
    for(i in 1:n)
    {
        f <- substr(tags[i], 1, 2)
        mtc <- match(f, tagLabs)
        if(is.na(mtc))
            stop("Field identifier '", f, "' not in List!\n")
        .Object@l[[f]] <- substr(tags[i], 4, nchar(tags[i]))
    }
    return(.Object)
})

setMethod("getHeaderText", signature="headerProgram", 
                                definition=function(object, delim="\t")
{
    n <- length(object@l)
    if(n==0)
        return(character(0))
    
    rfstr <- character(n)
    for(i in 1:n)
        rfstr[i] <- paste(names(object@l)[i], object@l[[i]], sep=":")
    return(paste("@PG", paste(rfstr, collapse=delim), sep=delim))
})

setMethod("getVal", signature="headerProgram", 
                                    definition=function(object, member)
{
    if(!is.character(member))
        stop("[getVal.headerProgram] Member must be character!\n")
    
    tagLabs <- c("ID", "PN", "CL", "PP", "DS", "VN")
    mtc <- match(member[1], tagLabs)
    
    if(is.na(mtc))
        stop("[getVal.headerProgram] Invalid member name!\n")
    
    return(object@l[[member]])
})

setMethod("setVal", signature="headerProgram", 
                                definition=function(object, members, values)
{
    if(!is.character(members) || !is.character(values))
        stop("Member name and value must be character!\n")
    
    if(length(members)!=length(values))
        stop("Members and values must have same length!\n")
    
    tagLabs <- c("ID", "PN", "CL", "PP", "DS", "VN")
    mtc <- match(members, tagLabs)
    if(any(is.na(mtc)))
        stop("Members must be valid Program Entries (See SAM Format Specification 1.3!\n")
    
    n <- length(members)
    obj <- deparse(substitute(object))
    for(i in 1:n)
    {
        txt <- paste(obj, "@l$", members[i], "<-'", values[i], "'", sep="")
        eval.parent(parse(text=txt))
    }
    return(invisible())
})

setMethod("as.list", signature="headerProgram", 
                                    definition=function(x, ...){return(x@l)})

setMethod("show", "headerProgram", function(object)
{
    n <- length(object@l)
    if(n>0)
    {
        cat("An object of class \"", class(object), "\"\n", sep="")
        for(i in 1:length(object@l))
        {
        cat(names(object@l)[i], ":", object@l[[i]], "\n")
        }
    }else{
        cat("An empty object of class \"", class(object), "\"\n", sep="")
    }
    return(invisible())
})

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#   End headerProgram
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#   bamHeaderText: Represents and manages textual version of bamHeader
#   See SAM Format Specification (v1.4-r985)
#
#   Contains header Segments :
#    head  = headerLine        : @HD Header Line
#    dict  = refSeqDict        : @SQ Reference Sequence dictionary
#    group = headerReadGroup   : @RG Read Group
#    prog  = headerProgram     : @PG Program
#
#    TODO:
#    com   = headerComment     : @CO One-line text comment
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#  Class definition and creational routines for bamHeaderText
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #


setMethod(f="initialize", signature="bamHeaderText", 
                                definition=function(.Object, bh="", delim="\n")
{
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    #  Parses Header data (as reported by getHeaderText)
    #  See Sam Format Specification 1.3 (Header Section)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    
    if(!is.character(bh))
        stop("[bamHeaderText.initialize] Argument must be string.\n")
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    #  Create empty header Set (so it's legal to call getHeaderText()')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    if(length(bh)==1 && nchar(bh)==0)
    {
        .Object@head <- new("headerLine")
        .Object@dict <- new("refSeqDict")
        .Object@group <- new("headerReadGroup")
        .Object@prog <- new("headerProgram")
        return(.Object)
    }
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    #  Split input string: Each fragment contains data for one header segment
    # 
    #  Identification of tags must be restricted on prefix 
    #  because there may be @RG entries present inside program segment
    #  (used as command line argument e.g. for aligner)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    bht <- unlist(strsplit(bh, split=delim))
    bht_pre <- substr(bht,1,3)
    
    # Read Header Line
    bhl <- bht[grep("@HD", bht_pre)]
    .Object@head <- new("headerLine", bhl)
    
    # Read Sequence Directory
    bsd <- bht[grep("@SQ", bht_pre)]
    .Object@dict <- new("refSeqDict", bsd)
    
    # Read Group
    brg <- bht[grep("@RG", bht_pre)]
    .Object@group <- new("headerReadGroup", brg)
    
    # Read Program Data
    bpd <- bht[grep("@PG", bht_pre)]
    .Object@prog <- new("headerProgram", bpd)
    
    # Read Text comment
    btc <- bht[grep("@CO", bht_pre)]
    com <- substring(btc, 3)
    return(.Object)
})

bamHeaderText <- function(head=NULL, dict=NULL, group=NULL, prog=NULL, com=NULL)
{
    bh <- new("bamHeaderText")
    if(!is.null(head))
    {
        if(is(head, "headerLine"))
            bh@head <- head
        else
            stop("[bamHeaderText] head must be 'headerLine'!")
    }
    if(!is.null(dict))
    {
        if(is(dict, "refSeqDict"))
            bh@dict <- dict
        else
        stop("[bamHeaderText] dict must be 'refSeqDict'")
    }
    
    if(!is.null(group))
    {
        if(is(group, "headerReadGroup"))
            bh@group <- group
        else
            stop("[bamHeaderText] group must be 'headerReadGroup'!")
    }
    
    if(!is.null(prog))
    {
        if(is(prog, "headerProgram"))
            bh@prog <- prog
        else
            stop("[bamHeaderText] prog must be 'headerProgram'!")
    }
    
    if(!is.null(com))
    {
        if(is.character(com))
            bh@com <- com
        else
            stop("[bamHeaderText] com must be 'character'!")
    }
    return(invisible(bh))
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#   End: Class definition and creational routines for bamHeaderText
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#   Public accessors for member objects for bamHeaderText
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
setMethod(f="headerLine", signature="bamHeaderText", 
                            definition=function(object) {return(object@head)})

setMethod(f="refSeqDict", signature="bamHeaderText", 
            definition=function(object) {return(object@dict)})



setMethod(f="headerReadGroup", signature="bamHeaderText", 
                            definition=function(object){return(object@group)})


setMethod(f="headerProgram", signature="bamHeaderText", 
                            definition=function(object){return(object@prog)})


setReplaceMethod("headerLine", "bamHeaderText", function(object, value)
{
    if(!is(value, "headerLine"))
        stop("[headerLine<-.bamHeaderText] value must be 'headerLine'!")
    object@head <- value
    return(object)
})


setReplaceMethod("refSeqDict", "bamHeaderText", function(object, value)
{
    if(!is(value, "refSeqDict"))
        stop("[refSeqDict<-.bamHeaderText] value must be 'refSeqDict'!")
    object@dict <- value
    return(object)
})


setReplaceMethod("headerReadGroup", "bamHeaderText", function(object, value)
{
    if(!is(value, "headerReadGroup"))
        stop("value must be 'headerReadGroup'!")
    object@group <- value
    return(object)
})


setReplaceMethod("headerProgram", "bamHeaderText", function(object, value)
{
    if(!is(value, "headerProgram"))
        stop("value must be 'headerProgram'!")
    object@prog <- value
    return(object)
})

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#   End: Public accessors for member objects for bamHeaderText
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #



setMethod("getHeaderText", signature="bamHeaderText",
                                    definition=function(object, delim="\n")
{
    hd <- getHeaderText(object@head)
    if(length(hd)==0)
        return(character(0))
    hd <- paste(hd, delim, sep="")
    
    dt <- getHeaderText(object@dict)
    if(length(dt)==0)
        return(character(0))
    dt <- paste(dt, delim, sep="")
    
    gp <- getHeaderText(object@group)
    if(length(gp)>0)
        gp <- paste(gp, delim, sep="")
    
    pg <- getHeaderText(object@prog)
    if(length(pg)>0)
        pg <- paste(pg, delim, sep="")
    
    if(length(object@com)>0)
        cm <- paste(paste("@CO", object@com, sep="\t"), collapse=delim)
    else
        cm <- character(0)
    return(paste(hd, dt, gp, pg, cm, sep=""))
})

setMethod("bamHeader", "bamHeaderText",  function(object)
{
    return(new("bamHeader", .Call("init_bam_header", getHeaderText(object))))
})


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#
#  bamWriter class
#  Encapsulates an write-opened Connection to a BAM-file.
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #


setMethod(f="initialize",  signature="bamWriter", 
          definition=function(.Object, header, filename)
{
    if(!is(header, "bamHeader"))
        stop("[initialize.bamWriter] header must be bamHeader!\n")
    
    if(!is.character(filename))
        stop("[initialize.bamWriter] filename must be character!\n")
    
    .Object@filename <- filename
    .Object@writer <- .Call("bam_writer_open", header@header, 
                                path.expand(filename), PACKAGE="rbamtools")
    
    return(.Object)
})

setMethod("filename",  "bamWriter",  function(object) return(object@filename))

setMethod("isOpen", signature="bamWriter", definition=function(con, rw="")
{
    return(!(.Call("is_nil_externalptr", con@writer, PACKAGE="rbamtools")))
})

setMethod(f="bamClose", signature="bamWriter", definition=function(object)
{
    return(invisible(.Call("bam_writer_close", 
                                        object@writer, PACKAGE="rbamtools")))
})

setMethod(f="bamSave", signature="bamWriter", 
                                definition=function(object, value, refid) 
{
    if(missing(refid))
        stop("[bamSave] refid is not optional!")
    
    if(!is.numeric(refid))
        stop("[bamSave] refid must be numeric")
    
    if(refid < 0)
        stop("[bamSave] refid must be >=0!")
    
    refid <- as.integer(refid)
    
    if(is(value, "bamAlign"))
    {
        return(invisible(.Call("bam_writer_save_align", object@writer, 
                            value@align, refid, PACKAGE="rbamtools")))
    }
    
    if(is(value, "bamRange")){
        return(invisible(.Call("bam_range_write", object@writer, 
                               value@range, refid, PACKAGE="rbamtools")))
    }
    else
        stop("bamSave: Saved object must be of type bamAlign or bamRange!\n")
})


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#
#  gapList
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

setMethod(f="initialize", "gapList", 
                definition=function(.Object, reader, coords, verbose=FALSE)
{
    if(!is(reader, "bamReader"))
    {
        cat("[initialize.gapList] Class of reader: ", class(reader), ".\n")
        stop("reader must be an instance of bamReader!\n")
    }
    
    if(length(coords) != 3)
        stop("coords must be 3-dim numeric (ref, start, stop)!\n")
    
    if(is.null(reader@index))
        stop("bamReader must have initialized index!\n")
    
    .Object@list <- .Call("gap_list_fetch", 
                reader@reader, reader@index, trunc(coords), PACKAGE="rbamtools")
    
    glsize <- .Call("gap_list_get_size", .Object@list, PACKAGE="rbamtools")
    
    if(verbose)
    {
        message("[initialize.gapList] Fetched list of size ", 
        format(glsize, big.mark=Sys.localeconv()[7]),
                                        " for refid ", coords[1], ".")
    }
    return(.Object)
})

# gapList function for retrieving objects in bamReader section

setMethod("size", signature="gapList", definition=function(object)
    {.Call("gap_list_get_size", object@list, PACKAGE="rbamtools")})

setMethod("nAligns", signature="gapList", definition=function(object)
    {.Call("gap_list_get_nAligns", object@list, PACKAGE="rbamtools")})

setMethod("nAlignGaps", signature="gapList", definition=function(object)
    {.Call("gap_list_get_nAlignGaps", object@list, PACKAGE="rbamtools")})

setMethod("show", "gapList", function(object)
{
    cat("An object of class '", class(object), "'. size: ",
                                                    size(object), "\n", sep="")
    
    cat("nAligns:", nAligns(object), "\tnAlignGaps:", nAlignGaps(object), "\n")
    return(invisible())
})

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#
#  gapSiteList
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #


setMethod(f="initialize", "gapSiteList", 
                                definition=function(.Object, reader, coords)
{
    if(missing(reader) || missing(coords))
        return(.Object)
    
    if(!is(reader, "bamReader"))
        stop("reader must be an instance of bamReader!\n")
    
    if(length(coords) != 3)
        stop("coords must be 3-dim numeric (ref, start, stop)!\n")
    
    if(is.null(reader@index))
        stop("bamReader must have initialized index!\n")
    
    .Object@list <- .Call("gap_site_list_fetch", 
        reader@reader, reader@index, trunc(coords), PACKAGE="rbamtools")
    
    return(.Object)
})

setMethod("size", signature="gapSiteList", definition=function(object)
{.Call("gap_site_list_get_size", object@list, PACKAGE="rbamtools")})

setMethod("nAligns", signature="gapSiteList", definition=function(object)
{.Call("gap_site_list_get_nAligns", object@list, PACKAGE="rbamtools")})

setMethod("nAlignGaps", signature="gapSiteList", definition=function(object)
{.Call("gap_site_list_get_nAlignGaps", object@list, PACKAGE="rbamtools")})

setMethod("refID", signature="gapSiteList", definition=function(object)
{.Call("gap_site_list_get_ref_id", object@list, PACKAGE="rbamtools")})


setMethod("show",  "gapSiteList", function(object)
{
    cat("An object of class '", class(object), 
                        "'. size: ", size(object), "\n", sep="")
    
    cat("nAligns:", nAligns(object), "\tnAlignGaps:", nAlignGaps(object), "\n")
    
    return(invisible())
})

merge.gapSiteList <- function(x, y, ...)
{
    if(!is(y,"gapSiteList"))
        stop("'y' must be of class 'gapSiteList'!")
    
    res <- new("gapSiteList")
    xref <- refID(x)
    
    if(refID(x) != refID(y))
        warning("[merge] 'x' and 'y' have different refID's. Using refID(x)!")
    
    res@list <- .Call("gap_site_list_merge", x@list, y@list, refID(x))
    
    return(res)
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#
#  bamGapList
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #


setMethod(f="initialize", "bamGapList", definition=function(.Object, reader)
{
    if(missing(reader))
    {
        .Object@list <- .Call("gap_site_ll_init")
        return(.Object)
    }
    
    if(!is(reader, "bamReader"))
        stop("reader must be an instance of bamReader!\n")
    if(is.null(reader@index))
        stop("bamReader must have initialized index!\n")
    
    ref <- getRefData(reader)
    ref$start <- 0L
    .Object@list <- .Call("gap_site_ll_fetch", 
                        reader@reader, reader@index, ref$ID, ref$start,
                        ref$LN, PACKAGE="rbamtools")
    
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    #  filter refdata for existing lists
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    sm <- .Call("gap_site_ll_get_summary_df", .Object@list, PACKAGE="rbamtools")
    mtc <- match(ref$ID, sm$ID)
    .Object@refdata <- ref[!is.na(mtc), ]
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    #  Re-enumerate ID's to 1:n
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    .Object@refdata$ID <- .Call("gap_site_ll_reset_refid", 
                                        .Object@list, PACKAGE="rbamtools")
    
    # ToDo: merge refdata with summary df?
  
    return(.Object)
})



setMethod("size", signature="bamGapList", definition=function(object)
    {.Call("gap_site_ll_get_size", object@list, PACKAGE="rbamtools")})

setMethod("nAligns", signature="bamGapList", definition=function(object)
    {.Call("gap_site_ll_get_nAligns", object@list, PACKAGE="rbamtools")})

setMethod("nAlignGaps", signature="bamGapList", definition=function(object)
    {.Call("gap_site_ll_get_nAlignGaps", object@list, PACKAGE="rbamtools")})

setMethod("show", "bamGapList", function(object)
{
    bm <- Sys.localeconv()[7]
    
    cat("An object of class '", class(object), "'. size: ", 
                    format(size(object), big.mark=bm), "\n", sep="")
    
    cat("nAligns:", format(nAligns(object), big.mark=bm), 
            "\tnAlignGaps:", format(nAlignGaps(object), big.mark=bm), "\n")
    
    return(invisible())
})

summary.bamGapList <- function(object, ...)
{
    return(merge(object@refdata, 
                        .Call("gap_site_ll_get_summary_df", object@list)))
}

merge.bamGapList <- function(x, y, ...)
{
    if(!is(y, "bamGapList"))
        stop("[merge.bamGapList] y must be bamGapList!")
    
    if(size(x)==0)
        stop("[merge.bamGapList] size(x)==0!")
    
    if(size(y)==0)
        stop("[merge.bamGapList] size(y)==0!")
    
    mref <- merge(x@refdata, y@refdata, by="SN", all=T)
    
    n <- dim(mref)[1]
    .Call("gap_site_ll_set_curr_first", x@list)
    .Call("gap_site_ll_set_curr_first", y@list)
    
    res <- new("bamGapList")
    for(i in 1:n)
    {
        if(is.na(mref$ID.x[i]))
        {
            .Call("gap_site_ll_add_curr_pp", 
                        y@list, res@list, as.integer(i-1))
            # copy values from .y to .x side (for later use in ref)
            mref[i, 2:4] <- mref[i, 5:7]
        }
        else if(is.na(mref$ID.y[i]))
        {
            .Call("gap_site_ll_add_curr_pp", 
                        x@list, res@list, as.integer(i-1))
        }else{
            .Call("gap_site_ll_add_merge_pp", 
                        x@list, y@list, res@list, as.integer(i-1))
        }
    }
    
    # get l-part of refdata
    ref <- mref[, 1:4]
    names(ref) <- c("SN", "ID", "LN", "start")
    
    # reset ID to new values  
    ref$ID <- 0:(n-1)
    res@refdata <- ref
    return(res)
}

readPooledBamGaps <- function(infiles, idxInfiles=paste(infiles, ".bai", sep=""))
{
    if(!is.character(infiles))
        stop("infiles must be character!")
    
    if(any(!file.exists(infiles)))
        stop("Files (infiles) not found!")
    
    if(!is.character(idxInfiles))
        stop("idxInFiles must be character!")
    
    if(any(!file.exists(idxInfiles)))
        stop("idxInfiles not found!")
    
    n <- length(infiles)  
    if(length(idxInfiles)!=n)
        stop("infiles and idxInfiles must have same length!")
    
    bm <- Sys.localeconv()[7]
    for(i in 1:n)
    {
        bam <- infiles[i]
        reader <- bamReader(bam)
        if(!file.exists(idxInfiles[i]))
        {
            message("[readPooledBamGaps] Creating BAM-index.", appendLF=FALSE)
            createIndex(reader, idxInfiles[i])
            message("Finished.")
        }
        
        loadIndex(reader, idxInfiles[i])
        message("[readPooledBamGaps] (",
                    format(i, width=2), "/", n, ")", appendLF=FALSE)
        
        if(i==1)
        ga <- bamGapList(reader)
        else
        {
            ga1 <- bamGapList(reader)
            ga <- merge.bamGapList(ga, ga1)
        }
            message("\tList-size: ", format(size(ga), width=7, big.mark=bm), 
                "\tnAligns: ", format(nAligns(ga), width=13, big.mark=bm), ".")
    }
    
    message("[readPooledBamGaps] Finished.")
    return(ga)
}

readPooledBamGapDf <- function(infiles, idxInfiles=paste(infiles, ".bai", sep=""))
{ 
    ga <- readPooledBamGaps(infiles, idxInfiles=paste(infiles, ".bai", sep=""))
    dfr <- as.data.frame(ga)
    attr(dfr, "nAligns") <- nAligns(ga)
    attr(dfr, "nAlignGaps") <- nAlignGaps(ga)
    return(dfr)
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#
#  bamRange
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#  Encapsulates a bunch of Alignment datasets that typically have been
#  read from a defined reference region in a BAM-file.
#  Technically,  the alignments are stored in a (C-implemented) double linked
#  list.
#  bamRange objects can be created by a reading procedure on an indexed
#  BAM-file. The alignments can be iterated, readed, written, deleted and
#  added. bamRange objects can be written to a BAM-file via an Instance
#  of bamWriter.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#  bamRange parameters:
#  1: seqid      : 0-based index of seqid
#  2: qrBegin    : 0-based left boundary of query region (query range begin)
#  3: qrEnd      : 0-based right boundary of query region (query range end)
#  4: complex    : 0= all aligns included, 1= only aligns with n_cigar > 1
#                                                    included
#  5: rSeqLen    : Length of reference sequence (from getRefData)
#  6: qSeqMinLen : Minimum of query sequence length (= read length)
#  7: qSeqMaxLen : Maximum of query sequence length (= read length)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #


bamRange <- function(object=NULL, coords=NULL, complex=FALSE)
{
    if(is.null(object))
        return(new("bamRange", NULL, NULL, FALSE))
    
    if(!is(object, "bamReader"))
        stop("object must be of class 'bamReader'!")
    
    if(!indexInitialized(object))
        stop("reader must have initialized index! Use 'loadIndex'!")

    return(new("bamRange", object, coords, complex))
}


setMethod(f="initialize", signature="bamRange", 
          definition=function(.Object, reader=NULL, coords=NULL, complex=FALSE)
{ 

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    #   Create empty range
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    if(is.null(reader))
    {
        .Object@range <- .Call("bam_range_init", PACKAGE="rbamtools")
        return(.Object)
    }
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    #   Create range from bam-file
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    if(!is(reader, "bamReader"))
        stop("reader must be an instance of bamReader!")
    
    # coords may either be missing or 3 entries are needed
    if(!is.null(coords))
    {
        if(length(coords)!=3)
            stop("coords must be numeric with length=3 (ref, start, stop)!")
    }

    
    if(is.null(reader@index))
        stop("bamReader must have initialized index!")
    
    if(!is(complex, "logical"))
        stop("complex must be logical!")
    
    if(length(complex) > 1)
        stop("complex must have length 1!")
    
    if(!indexInitialized(reader))
        stop("reader must have initialized index! Use 'loadIndex'!")
    
    .Object <- .Call("bam_range_fetch", reader@reader, 
                    reader@index, trunc(coords), complex, PACKAGE="rbamtools")
    
    return(.Object)
})


setMethod("size", signature="bamRange", definition=function(object)
    {.Call("bam_range_get_size", object@range, PACKAGE="rbamtools")})


setMethod("getCoords", "bamRange", function(object)
    { return(.Call("bam_range_get_coords", object@range))})

setMethod("getParams", "bamRange", function(object)
    { return(.Call("bam_range_get_params", object@range))})

setMethod("getSeqLen", "bamRange", function(object){
  return(.Call("bam_range_get_seqlen", object@range, PACKAGE="rbamtools"))
})


setMethod("getRefName", "bamRange", function(object)
    return(.Call("bam_range_get_refname", object@range, PACKAGE="rbamtools")))


setMethod("show", "bamRange", function(object){
  bm <- Sys.localeconv()[7]
  w <- 11
  r <- "right"
  cat("Class       : ", format(class(object), w=w, j=r)                   , "\n", sep="")
  cat("Size        : ", format(format(size(object), big.m=bm), w=w, j=r)   , "\n", sep="")

  params <- .Call("bam_range_get_params", object@range, PACKAGE="rbamtools")
  cat("Seqid       : ", format(format(params[1], big.m=bm), w=w, j=r)     , "\n", sep="")
  cat("qrBegin     : ", format(format(params[2], big.m=bm), w=w, j=r)     , "\n", sep="")
  cat("qrEnd       : ", format(format(params[3], big.m=bm), w=w, j=r)     , "\n", sep="")
  cat("Complex     : ", format(params[4], w=w, big.m=bm)                , "\n", sep="")
  cat("rSeqLen(LN) : ", format(format(params[5], big.m=bm), w=w, j=r)   , "\n", sep="")
  cat("qSeqMinLen  : ", format(format(params[6], big.m=bm), w=w, j=r)   , "\n", sep="")
  cat("qSeqMaxLen  : ", format(format(params[7], big.m=bm), w=w, j=r)   , "\n", sep="")
  
  refname <- .Call("bam_range_get_refname", object@range, PACKAGE="rbamtools")
  if(!is.null(refname))
  cat("Refname     : ", format(refname, w=w, j="right")       , "\n", sep="")  
  return(invisible())
})


setMethod("getAlignRange", "bamRange", function(object)
    return(.Call("bam_range_get_align_range", 
                object@range, PACKAGE="rbamtools")))


setMethod("getNextAlign", signature="bamRange", definition=function(object)
{
    ans <- .Call("bam_range_get_next_align", object@range, PACKAGE="rbamtools")
    # Must be checked because align list returns NULL when end is reached
    if(is.null(ans))
        return(ans)
    else
        return(new("bamAlign", ans))
})


setMethod("getPrevAlign", signature="bamRange", definition=function(object)
{
    return(new("bamAlign", .Call("bam_range_get_prev_align", 
                        object@range, PACKAGE="rbamtools")))
})


setMethod("stepNextAlign", signature("bamRange"), definition=function(object)
{
    .Call("bam_range_step_next_align", object@range)
    return(invisible())
})


setMethod("stepPrevAlign", signature("bamRange"), definition=function(object)
{
    .Call("bam_range_step_prev_align", object@range)
    return(invisible())
})


# Resets current align to NULL position (i.e. before first element)
# The next call to getNextAlign then returns the first element of list
setMethod("rewind", signature="bamRange", definition=function(object)
{
    invisible(.Call("bam_range_wind_back", object@range, PACKAGE="rbamtools"))
})


setMethod("push_back", signature="bamRange", definition=function(object, value)
{
    if(!is(value, "bamAlign"))
        stop("pushed object must be of class \"bamAlign\"\n")
    
    .Call("bam_range_push_back", object@range, value@align, PACKAGE="rbamtools")
})


setMethod("pop_back", signature="bamRange", definition=function(object)
{
    .Call("bam_range_pop_back", object@range, PACKAGE="rbamtools")
})


setMethod("push_front", signature="bamRange", definition=function(object, value)
{
    if(!is(value, "bamAlign"))
        stop("pushed object must be of class \"bamAlign\"\n")
    
    .Call("bam_range_push_front", object@range, value@align, PACKAGE="rbamtools")
})


setMethod("pop_front", signature="bamRange", definition=function(object)
{
    .Call("bam_range_pop_front", object@range, PACKAGE="rbamtools")
})


setMethod("writeCurrentAlign", signature="bamRange", definition=function(object, value)
{
    if(!is(value, "bamAlign"))
        stop("written object must be of class \"bamAlign\"\n")
    
    .Call("bam_range_write_current_align",
                        object@range, value@align, PACKAGE="rbamtools")
})


setMethod("insertPastCurrent", signature="bamRange", 
                                            definition=function(object, value)
{
    if(!is(value, "bamAlign"))
        stop("written object must be of class \"bamAlign\"\n")
    
    .Call("bam_range_insert_past_curr_align", 
                                object@range, value@align, PACKAGE="rbamtools")
})


setMethod("insertPreCurrent", signature="bamRange",
                                            definition=function(object, value)
{
    if(!is(value, "bamAlign"))
        stop("written object must be of class \"bamAlign\"\n")
    
    .Call("bam_range_insert_pre_curr_align", 
                        object@range, value@align, PACKAGE="rbamtools")
})


setMethod("moveCurrentAlign", signature="bamRange", 
                                        definition=function(object, target)
{
    if(!is(target, "bamRange"))
        stop("target must be bamRange!\n")
    
    .Call("bam_range_mv_curr_align", object@range, target@range)
    return(invisible())
})


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#  Rudimentary subsetting operator:
#  Does not change the order of elements, just returns subset
#  Therefore sorts given index i
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

setMethod("[", signature="bamRange", function(x, i)
{
    i <- sort(as.integer(i))
    if(i[1] < 1)
        stop("No negative indices allowed. Use 'pop_front' or 'pop_back'!")
    
    if(i[length(i)] > size(x))
        stop("Out of bounds index (> size(x))!")
    
    return(.Call("bam_range_idx_copy", x@range, i, PACKAGE="rbamtools"))
})


setMethod("range2fastq", signature="bamRange",
    definition=function(object, filename, which, append=FALSE)
{
    message("[range2fastq] Function is deprecated. Use rangeToFastq.")
    return(rangeToFastq(object, filename, which, append))
})

setMethod("rangeToFastq", signature="bamRange",
    definition=function(object, filename, which, append=FALSE)
{
    if(!is.character(filename))
        stop("'filename' must be character!")
    
    if(!is.logical(append))
        stop("'append' must be logical!")
    
    if(missing(which))
    {
        .Call("bam_range_write_fastq", object@range, filename,
                            append, PACKAGE="rbamtools")
    }else{
        if(!is.numeric(which))
            stop("'which' must be numeric!")
        
        mx <- max(which)
        
        if(mx>size(object))
        {
            cat("[rangeToFastq] Maximum index (", mx,
                    ") is greater than size of range (",
                    size(object), ")!\n", sep="")
        }
        
        written <- .Call("bam_range_write_fastq_index",
            object@range, filename, as.integer(sort(unique(which))),
            append, PACKAGE="rbamtools")
        
        cat("[rangeToFastq]", written, "records written.\n")
    }
    return(invisible())
})

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#  Functions to read and display phred qualities from bamRange
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

setMethod("getQualDf", "bamRange", function(object, prob=FALSE, ...)
{
    if(!is.logical(prob))
        stop("[getQualDf.bamRange] ")
    if(prob)
    {
        qdf <- .Call("bam_range_get_qual_df", 
                                        object@range, PACKAGE="rbamtools")
        
        rel <- function(x)
        {
            xs <- sum(x)
            if(xs > 0)
                return(x / xs)
            return(x)
        }
        
        res <- data.frame(lapply(qdf, rel))
        names(res) <- names(qdf)
        attributes(res)$col.sums <- unlist(lapply(qdf, sum))
        return(res)
    }
    return(.Call("bam_range_get_qual_df", object@range, PACKAGE="rbamtools"))
})


setMethod("getQualQuantiles", "bamRange", function(object, quantiles, ...)
{
    if(!is.numeric(quantiles))
        stop("[getQualQuantiles.bamRange] quantiles must be numeric!")
    
    if(!(all(quantiles >= 0) & all(quantiles <= 1)))
        stop("[getQualQuantiles.bamRange] all quantiles mustbe in [0, 1]")
    
    quantiles <- sort(unique(round(quantiles, 2)))

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    #  Count qual values for each sequence position
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    qdf <- .Call("bam_range_get_qual_df", object@range, PACKAGE="rbamtools")
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    #  Convert integer counts into column-wise relative values
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    rel <- function(x)
    {
        xs <- sum(x)
        if(xs>0)
            return(x/xs)
        return(x)
    }
    
    qrel <- data.frame(lapply(qdf, rel))
    names(qrel) <- names(qdf)
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    #  Walk through each column and extract row number
    #  for given quantile values
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    res <- .Call("get_col_quantiles", quantiles, qrel, PACKAGE="rbamtools")
    return(res)
})


setMethod("plotQualQuant",  "bamRange",  function(object)
{
    quant <- c(0.1, 0.25, 0.5, 0.75, 0.9)
    cols <- c("#1F78B4", "#FF7F00", "#E31A1C", "#FF7F00", "#1F78B4")
    
    qq <- getQualQuantiles(object, quant)
    
    maxQ <- floor(1.2 * max(qq))
    xv <- 1:ncol(qq)
    
    plot(xv, xv, ylim=c(0, maxQ), type="n", bty="n", las=1,
            ylab="phred score", xlab="sequence position",
            main="Phred Quantiles for sequence")
    
    lines(xv, qq[1, ], col=cols[1], lty=2)
    lines(xv, qq[2, ], col=cols[2], lty=1)
    lines(xv, qq[3, ], col=cols[3], lwd=2)
    lines(xv, qq[4, ], col=cols[4], lty=1)
    lines(xv, qq[5, ], col=cols[5], lty=2)
    
    legend("top", ncol=6, lty=c(2, 1, 1, 1, 2),
        lwd=c(1, 1, 2, 1, 1), col=cols, xjust=0.5,
        legend=c("10%", "25%", "50%", "75%", "90%"), bty="n", cex=0.8)
    
    return(invisible()) 
})



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#
#  alignDepth
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#  alignDepth parameters:
#  - - bamRange derived - -
#  1: seqid      : 0-based index of seqid
#  2: qrBegin    : 0-based left boundary of query region (query range begin)
#  3: qrEnd      : 0-based right boundary of query region (query range end)
#  4: complex    : 0= all aligns included, 1= only aligns with n_cigar > 1
#                                                    included
#  5: rSeqLen    : Length of reference sequence (from getRefData)
#  6: qSeqMinLen : Minimum of query sequence length (= read length)
#  7: qSeqMaxLen : Maximum of query sequence length (= read length)
#  - - alignDepth proprietary - -
#  6: gap     : 0=all aligns counted, 1=only gap adjacent match regions
#  counted
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #



setMethod("alignDepth", "bamRange", function(object, gap=FALSE)
{
    if(!is.logical(gap))
        stop("gap must be logical!")
    
    return(.Call("bam_range_get_align_depth", 
                        object@range, gap, PACKAGE="rbamtools"))
})

setMethod("show", "alignDepth", function(object){
    bm <- Sys.localeconv()[7]
    w <- 11
    cat("Class       : ", format(class(object), w=w, j="right")  , "\n", sep="")
    cat("Seqid       : ", format(object@params[1], w=w, big.m=bm)    , "\n", sep="")
    cat("qrBegin     : ", format(object@params[2], w=w, big.m=bm)    , "\n", sep="")
    cat("qrEnd       : ", format(object@params[3], w=w, big.m=bm)    , "\n", sep="")
    cat("Complex     : ", format(object@params[4], w=w, big.m=bm)    , "\n", sep="")
    cat("rSeqLen(LN) : ", format(object@params[5], w=w, big.m=bm)    , "\n", sep="")
    cat("qSeqMinLen  : ", format(object@params[6], w=w, big.m=bm)    , "\n", sep="")  
    cat("qSeqMaxLen  : ", format(object@params[7], w=w, big.m=bm)    , "\n", sep="")
    cat("refname     : ", format(object@refname, w=w, j="right")     , "\n", sep="") 
    n <- 6
    x <- object@depth[1:n]
    names(x) <- object@pos[1:n]
    print(x)
    return(invisible())
})


setMethod("getDepth", "alignDepth", function(object, named=FALSE)
{
    if(!is.logical(named))
        stop("[getDepth.alignDepth] named must be logical!")
    if(named)
    {
        dp <- object@depth
        names(dp)=object@pos
        return(dp)
    }
    return(object@depth)
})


setMethod("getPos",     "alignDepth", function(object) {return(object@pos)})
setMethod("getParams", "alignDepth", function(object) {return(object@params)})


setMethod("plotAlignDepth", "alignDepth",
        function(object, start=NULL, end=NULL, xlim=NULL,
                    main="Align Depth", xlab="Position", 
                    ylab="Align Depth",  transcript="",
                    strand=NULL , log="y", cex.main=2,
                    col="grey50", fill="grey90", grid=TRUE, 
                    box.col="grey20", box.border="grey80", ... )
{
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    #  Start and end positions for exon - rectangles
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    
    if(!is.null(start))
    {
        if(is.null(end))
            stop("'end' must be given when start is present")
        
        if(length(start) != length(end))
            stop("'start' and 'end' must have same length!")
        
        if(any(start <0 ) || any(end < 0))
           stop("No negative values allowed in 'start' and 'end'")
    }
    
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    #  Prepare align depth values for polygon plotting and
    #  logarithmic scale
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    nx <- length(object@pos)
    x <- c(object@pos[1], object@pos, object@pos[nx])
    y <- c(1, ifelse(object@depth == 0, 1, object@depth), 1)
    
    
    #  
    if(is.null(xlim))
        xlim <- c(x[1], x[nx])
    
    
    if(!is.null(start))
    {
        #  + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + # 
        #  Prepare plot area and layout
        #  + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + # 
        
        # Reset values
        op <- par(no.readonly=TRUE)
        par(oma=c(3, 2, 1, 3))
        
        m <- matrix(c(1, 0, 2), ncol=1)
        layout(m, heights=c(5, lcm(0.2), 2))
        
        #  + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + # 
        #  Do upper plot: align depth polygon
        #  + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + # 
        par(mar=c(0, 4, 3, 1) + 0.1)
        
        plot(x, y, type="n", las=1, main=main, xlim=xlim, 
            xlab="", ylab="", log=log, xaxt="n", bty="n",
            cex.axis=1.5, cex.main=cex.main, ...)
        
        polygon(x, y, col=fill, border=col, lwd=2)
        
        if(grid)
            grid()
        
        #
        mtext(paste("Refname:", object@refname), adj=1, cex=0.8)
        # ylab has to be positioned more outside
        mtext(ylab, side=2, line=4 , adj=0.5)
        
        #  + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + # 
        #  Do lower plot: Draw exon boxes, x-axis and transcript text
        #  + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + # 
        par(mar=c(4, 4, 0, 1) + 0.1)
        plot(x, y, type="n", yaxt="n", xlim=xlim, ylim=c(0, 10), 
            cex.axis=1.5, bty="n", ylab="", xlab=xlab , 
            cex.lab=1.5, ...)
        
        # Draw horizontal gene - line
        if(is.null(strand))
            lines(xlim,c(8,8))
        else
        {
            if(strand == "+")
            {
                arrows(x0=xlim[1], y0=8, x1=xlim[2], y1=8, code=2)
            } else {
                arrows(x0=xlim[1], y0=8, x1=xlim[2], y1=8, code=1)
            }
        }
        
        
        # Draw exon boxes
        for(i in 1:length(start))
            rect(start[i], 6, end[i], 10, col=box.col, border=box.border)
        
        # Write transcript text
        text(xlim[1], 2, transcript, adj=0, cex=1.2)
        
        #  + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + # 
        #  Cleanup
        #  + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + # 
        par(op)
        
    }else{
        
        #  + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + # 
        #  Alternative plot when no exon positions are given
        #  + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + # 
        plot(x, y, type="l", las=1, col=col, bty="n" , log=log,
                    xlab=xlab, ylab=ylab, main=main, ...)
        # log="" turns log scaling off.
        
        if(grid)
            grid()
        mtext(paste("Refname:", object@refname))
    }
    
    return(invisible())
})


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#  Count nucleotides
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
setMethod("countNucs", "bamRange", function(object)
                { return(.Call("bam_range_count_nucs",
                            object@range, PACKAGE="rbamtools"))})


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#
#  bamAlign
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#  bamAlign encapsulates all contained data in a single dataset in
#  a BAM-file.
#  bamAlign objects can be read from a bamReader instance and written to a
#  bamWriter instance. All contained data can be read and written via
#  accessors functions.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #


setMethod(f="initialize", signature="bamAlign",
                                    definition=function(.Object,align=NULL)
{
    .Object@align <- align
    return(.Object)
})


setMethod("show", "bamAlign", function(object)
{
    bm <- Sys.localeconv()[7]
    w <- 11
    r <- "right"
    cat("Class       : ", format(class(object)  , w=w, j=r)                     , "\n", sep="")
    cat("refId       : ", format(refID(object)  , w=w, j=r)                     , "\n", sep="") 
    cat("Position    : ", format(format(position(object), big.m=bm), w=w, j=r)   , "\n", sep="")
    
    cat("\nCigar Data  :\n")
    print(cigarData(object))
})


bamAlign <- function(qname, qseq, qqual, cigar, refid, position, 
                flag=272L, alqual=10L, mrefid=(-1L), mpos=(-1L), insertsize=0L)
{
    if(missing(qname))
        stop("[bamAlign] Missing query name!")
    
    if(missing(qseq))
        stop("[bamAlign] Missing query sequence string!")
    
    if(missing(qqual))
        stop("[bamAlign] Missing query quality string!")
    
    if(missing(cigar))
        stop("[bamAlign] Missing CIGAR string!")
    
    if(missing(refid))
        stop("[bamAlign] Missing refid!")
    
    if(missing(position))
        stop("[bamAlign] Missing position!")
    
    
    if(!is.character(qname))
        stop("[bamAlign] Query name must be character!")
    
    if(!is.character(qseq))
        stop("[bamAlign] Query sequence must be character!")
    
    if(!is.character(qqual))
        stop("[bamAlign] Query quality must be character!")
    
    if(nchar(qseq)!=nchar(qqual))
        stop("Query sequence string and quality string must have equal size!")
    
    
    if(!is.character(cigar))
        stop("[bamAlign] CIGAR string must be character!")
    
    if(!is.numeric(refid))
        stop("[bamAlign] refid must be numeric!")
    
    if(!is.numeric(position))
        stop("[bamAlign] position must be numeric!")
    
    refid <- as.integer(refid)
    position <- as.integer(position)
    
    # String-values:
    # 1) query-name
    # 2) query sequence
    # 3) quality string
    # 4) CIGAR string
    
    strval <- character(4)
    strval[1] <- qname
    strval[2] <- qseq
    strval[3] <- qqual
    strval[4] <- cigar
    
    # Integer-values:
    # 1) refid
    # 2) position
    # 3) flag
    # 4) align quality
    # 5) mate refid
    # 6) mate position
    # 7) insert size
    intval <- integer(4)
    intval[1] <- refid
    intval[2] <- position
    intval[3] <- flag
    intval[4] <- alqual
    intval[5] <- mrefid
    intval[6] <- mpos
    intval[7] <- insertsize
    
    ans <- .Call("bam_align_create", strval, intval)
    
    # Must be checked because align list returns NULL when end is reached
    if(is.null(ans))
    {  
        cat("[bamAlign] Align creation unsuccessful! Data inconsistency?\n")
        return(NULL)
    }
    
    return(new("bamAlign", ans))
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#  bamAlign Member Reader functions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

setMethod(f="name", signature="bamAlign", definition=function(object)
{
    .Call("bam_align_get_name", object@align, PACKAGE="rbamtools")
})

setMethod(f="refID", signature="bamAlign", definition=function(object)
{
    .Call("bam_align_get_refid", object@align, PACKAGE="rbamtools")
})

setMethod(f="position", signature="bamAlign", definition=function(object)
{
    .Call("bam_align_get_position", object@align, PACKAGE="rbamtools")
})

setMethod("nCigar", signature="bamAlign", definition=function(object)
{
    .Call("bam_align_get_nCigar", object@align, PACKAGE="rbamtools")
})

setMethod(f="cigarData", signature="bamAlign", definition=function(object)
{
    .Call("bam_align_get_cigar_df", object@align, PACKAGE="rbamtools")
})

setMethod(f="mateRefID", signature="bamAlign", definition=function(object)
{
    .Call("bam_align_get_mate_refid", object@align, PACKAGE="rbamtools")
})

setMethod(f="matePosition", signature="bamAlign", definition=function(object)
{
    .Call("bam_align_get_mate_position", object@align, PACKAGE="rbamtools")
})

setMethod(f="insertSize", signature="bamAlign", definition=function(object)
{
    .Call("bam_align_get_insert_size", object@align, PACKAGE="rbamtools")
})

setMethod(f="mapQuality", signature="bamAlign", definition=function(object)
{
    .Call("bam_align_get_map_quality", object@align, PACKAGE="rbamtools")
})

setMethod(f="alignSeq", signature="bamAlign", definition=function(object)
{
    .Call("bam_align_get_segment_sequence", object@align, PACKAGE="rbamtools")
})

setMethod(f="alignQual", signature="bamAlign", definition=function(object)
{
    .Call("bam_align_get_qualities", object@align, PACKAGE="rbamtools")
})

setMethod(f="alignQualVal", signature="bamAlign", definition=function(object)
{
    .Call("bam_align_get_qual_values", object@align, PACKAGE="rbamtools")
})




# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#  1.4  The alignment section: mandatory fields
#       2. FLAG: bitwise FLAG
#  Queries against alignment flag (Readers and Accessors)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#  0x1 template having multiple segments in sequencing
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
setMethod("paired", "bamAlign", function(object)
{
    .Call("bam_align_is_paired", object@align, PACKAGE="rbamtools")
})

setReplaceMethod(f="paired",
    signature="bamAlign", definition=function(object, value)
{
    if(!is.logical(value))
        stop("value must be boolean")
    
    .Call("bam_align_set_is_paired", 
                    object@align, value, PACKAGE="rbamtools")
    return(object)
})


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#  0x2 each segment properly aligned according to the aligner
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
setMethod("properPair", "bamAlign", function(object)
{
    .Call("bam_align_mapped_in_proper_pair", object@align, PACKAGE="rbamtools")
})

setReplaceMethod(f="properPair",
    signature="bamAlign", definition=function(object, value)
{
    if(!is.logical(value))
        stop("value must be boolean")
    
    .Call("bam_align_set_mapped_in_proper_pair",
                    object@align, value, PACKAGE="rbamtools")
    
    return(object)
})


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#  0x4 segment unmapped
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
setMethod("unmapped", "bamAlign", function(object)
{
    .Call("bam_align_is_unmapped", object@align, PACKAGE="rbamtools")
})

setReplaceMethod(f="unmapped",
    signature="bamAlign", definition=function(object, value)
{
    if(!is.logical(value))
        stop("value must be boolean")
    
    .Call("bam_align_set_is_unmapped",
                    object@align, value, PACKAGE="rbamtools")
    return(object)
})


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#  0x8 next segment in the template unmapped
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
setMethod("mateUnmapped", "bamAlign", function(object)
{
    .Call("bam_align_mate_is_unmapped", object@align, PACKAGE="rbamtools")
})

setReplaceMethod(f="mateUnmapped",
    signature="bamAlign", definition=function(object, value)
{
    if(!is.logical(value))
        stop("value must be boolean")
    
    .Call("bam_align_set_mate_is_unmapped",
        object@align, value, PACKAGE="rbamtools")
    return(object)
})


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#  0x10 SEQ being reverse complemented
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
setMethod("reverseStrand", "bamAlign", function(object)
{
    .Call("bam_align_strand_reverse", object@align, PACKAGE="rbamtools")
})

setReplaceMethod(f="reverseStrand",
    signature="bamAlign", definition=function(object, value)
{
    if(!is.logical(value))
        stop("value must be boolean")
    
    .Call("bam_align_set_strand_reverse", 
        object@align, value, PACKAGE="rbamtools")
    return(object)
})


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#  0x20 SEQ of the next segment in the template being reversed
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
setMethod("mateReverseStrand", "bamAlign", function(object)
{
    .Call("bam_align_mate_strand_reverse", object@align, PACKAGE="rbamtools")
})

setReplaceMethod(f="mateReverseStrand",
    signature="bamAlign", definition=function(object, value)
{
    if(!is.logical(value))
        stop("value must be boolean")
    
    .Call("bam_align_set_mate_strand_reverse",
        object@align, value, PACKAGE="rbamtools")
    return(object)
})


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#  0x40 the first segment in the template
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
setMethod("firstInPair",  "bamAlign",  function(object)
{
    .Call("bam_align_is_first_in_pair", object@align, PACKAGE="rbamtools")
})

setReplaceMethod(f="firstInPair", 
    signature="bamAlign", definition=function(object, value)
{
    if(!is.logical(value))
        stop("class bamReader, FirstInPair setter: value must be boolean")
     
    .Call("bam_align_set_is_first_in_pair", 
            object@align, value, PACKAGE="rbamtools")
    return(object)
})


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#  0x80 the last segment in the template
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

setMethod("secondInPair", "bamAlign", function(object)
{
    .Call("bam_align_is_second_in_pair", object@align, PACKAGE="rbamtools")
})

setReplaceMethod(f="secondInPair",
    signature="bamAlign", definition=function(object, value)
{
    if(!is.logical(value))
        stop("value must be boolean")
    
    .Call("bam_align_set_is_second_in_pair", 
        object@align, value, PACKAGE="rbamtools")
                     
    return(object)
})


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#  0x100 secondary alignment
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
setMethod("secondaryAlign", "bamAlign", function(object)
{
  .Call("bam_align_is_secondary_align",object@align, PACKAGE="rbamtools")
})

setReplaceMethod(f="secondaryAlign", signature="bamAlign",
                                        definition=function(object, value)
{
    if(!is.logical(value))
        stop("value must be boolean")
    
    .Call("bam_align_set_is_secondary_align", 
                    object@align, value, PACKAGE="rbamtools")
    return(object)
})


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#  0x200 not passing quality controls
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
setMethod("failedQC", "bamAlign", function(object)
{
    return(.Call("bam_align_fail_qc", object@align, PACKAGE="rbamtools"))
})

setReplaceMethod(f="failedQC", 
        signature="bamAlign", definition=function(object, value)
{
    if(!is.logical(value))
        stop("class bamReader, failedQC setter: value must be boolean")
    
    .Call("bam_align_set_fail_qc", object@align, value, PACKAGE="rbamtools")
    return(object)
})


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#  0x400 PCR or optical duplicate
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
setMethod("pcrORopt_duplicate", "bamAlign", function(object)
{
    return(.Call("bam_align_is_pcr_or_optical_dup",
                object@align, PACKAGE="rbamtools"))
})

setReplaceMethod(f="pcrORopt_duplicate",
        signature="bamAlign", definition=function(object, value)
{
    if(!is.logical(value))
        stop("class bamReader, Duplicate setter: value must be boolean")
    .Call("bam_align_set_is_pcr_or_optical_dup", 
                            object@align, value, PACKAGE="rbamtools")
    return(object)
})


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#  0x800 supplementary alignment
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
setMethod("suppAlign", "bamAlign", function(object)
{
    .Call("bam_align_is_supplementary_align",object@align, PACKAGE="rbamtools")
})

setReplaceMethod(f="suppAlign",
        signature="bamAlign", definition=function(object, value)
{
    if(!is.logical(value))
        stop("value must be boolean")
    
    .Call("bam_align_set_is_supplementary_align", 
        object@align, value, PACKAGE="rbamtools")
    return(object)
})


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#  flag
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
setMethod("flag", "bamAlign", function(object)
{
  .Call("bam_align_get_flag", object@align, PACKAGE="rbamtools")
})

setReplaceMethod(f="flag", signature="bamAlign",
                                    definition=function(object, value)
{
    if(!is.numeric(value))
        stop("value must be numeric")
    
    if(!is.integer(value))
    {
        value <- as.integer(value)
        message("[flag] Value is coerced to integer.")
    }
    
    .Call("bam_align_set_flag", object@align, value, PACKAGE="rbamtools")
    return(object)
})


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#   End: Queries against alignment flag (Readers and Accessors)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
setMethod("countNucs","bamAlign",function(object)
{
    return(.Call("bam_align_count_nucs", object@align, PACKAGE="rbamtools"))
})


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#   End: bamAlign
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#   coercing
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

as.data.frame.bamRange <- function(x, row.names=NULL, optional=FALSE, ...)
{
    return(.Call("bam_range_get_align_df", x@range, PACKAGE="rbamtools"))
}

as.data.frame.gapList <- function(x, row.names=NULL, optional=FALSE, ...)
{
    return(.Call("gap_list_get_df", x@list, PACKAGE="rbamtools"))
}

as.data.frame.gapSiteList <- function(x, row.names=NULL, optional=FALSE, ...)
{
    return(.Call("gap_site_list_get_df", x@list, PACKAGE="rbamtools"))
}

as.data.frame.bamGapList <- function(x, row.names=NULL, optional=FALSE, ...)
{
    return(.Call("gap_site_ll_get_df", x@list,
                                        x@refdata$SN, PACKAGE="rbamtools"))
}

as.data.frame.refSeqDict <- function(x, row.names=NULL, optional=FALSE, ...)
{
    n <- length(x@SN)
    if(n == 0)
    {
        return(data.frame(SN=character(0), LN=numeric(0), AS=character(0),
                        M5=numeric(0), SP=character(0) , UR=character(0)))
    }
    
    if(is.null(row.names))
        row.names <- 1:(length(x@SN))
    else if(length(row.names) != length(x@SN))
        stop("length(row.names)!=length(x@SN)!")
    
    return(data.frame(SN=x@SN, LN=x@LN, 
                AS=x@AS ,M5=x@M5, SP=x@SP, UR=x@UR, row.names=row.names))
}


setAs("bamRange","data.frame", function(from)
{
    return(.Call("bam_range_get_align_df", from@range, PACKAGE="rbamtools"))
})

setAs("gapList", "data.frame", function(from)
{
    return(.Call("gap_list_get_df", from@list, PACKAGE="rbamtools"))
})

setAs("refSeqDict", "data.frame", function(from)
{
    return(data.frame(SN=from@SN, LN=from@LN, AS=from@AS, M5=from@M5,
                    SP=from@SP, UR=from@UR, row.names=1:length(from@SN)))
})

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#  End: coercing
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#  Miscellaneous functions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #


createIdxBatch <- function(bam, idx=paste(bam, ".bai", sep=""), rebuild=FALSE)
{
    if(!is.character(bam))
        stop("'bam' must be character!")
    
    if(!is.character(idx))
        stop("'idx' must be character!")
    
    if(length(bam)!=length(idx))
        stop("'bam' and 'idx' must have same length!")
    
    if(!is.logical(rebuild))
        stop("'rebuild' must be logical!")
    
    if(length(rebuild) > 1)
        stop("'rebuild' must have length 1!")
    
    
    n <- length(bam)
    for(i in 1:n)
    {
        message("[", format(i, width=2), "/", n, "] ", appendLF=FALSE)
        if(!file.exists(bam[i]))
        stop("File ", i, " does not exist!")
        
        if(rebuild[1])
        {
        reader <- bamReader(bam[i])
        createIndex(reader, idx[i])
        bamClose(reader)       
        }else{
            if(!file.exists(idx[i]))
            {
                reader <- bamReader(bam[i])
                createIndex(reader, idx[i])
                bamClose(reader)     
            }
        }
        message("OK.")
    }
    return(invisible())
}

create.idx.batch <- function(bam, idx=paste(bam, ".bai", sep=""), rebuild=FALSE)
{ .Deprecated("createIdxBatch",package="rbamtools") }



countTextLines <- function(filenames)
{
    if(!is.character(filenames))
        stop("[countTextLines] filename must be character!")
    
    if(!all(file.exists(filenames)))
        stop("[countTextLines] Missing files!")
    
    n <- length(filenames)
    res <- numeric(n)
    bm <- Sys.localeconv()[7]
    
    for(i in 1:n)
    {
        cat("[countTextLines] Counting '", basename(filenames[i]), "'", sep="")
        res[i] <- .Call("count_text_lines", filenames[i])
        cat("\t found:", format(res[i], big.mark=bm), ".\n")
    }
    cat("[countTextLines] Finished. Found", 
                                format(sum(res), big.mark=bm), "lines.\n")
    return(res)
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#   Unexported and undocumented routines
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

readSepGapTables <- function(bam, profo, defo="sep_gap",
                                            idx=paste(bam, ".bai", sep=""))
{
    require(rbamtools)
    fo <- file.path(profo, defo)
    
    if(!file.exists(fo))
        dir.create(fo)
    bm <- Sys.localeconv()[7]
    
    n <- length(bam)
    for(i in 1:n)
    {
        cat("[readSepGapTables] i:(", format(i, width=2), "/", n, ")", sep="")
        
        if(!file.exists(bam[i]))
            stop("[readSepGapTables] i:", i, " File does not exist!")
        
        reader <- bamReader(bam[i])
        
        if(!file.exists(idx[i]))
            createIndex(reader, idx[i])
        
        loadIndex(reader, idx[i])
        
        if(i==1)
        {
            bsl <- bamGapList(reader)
            dfr <- as.data.frame(bsl)
            save(dfr, file=file.path(fo, paste("bsl_", i, ".RData", sep="")))
            write.table(dfr, file=file.path(fo, 
                paste("bsl_", i, ".csv", sep="")), sep=";", row.names=FALSE)
            
            cat("\r[readSepGapTables] i:(", format(i, width=2),
                "/", n, ")\tnr sites: ", format(size(bsl), big.mark=bm,
                                                width=9), "\n", sep="")
            
        }else{
            # save site-table for bam[i]
            bsli <- bamGapList(reader)
            dfri <- as.data.frame(bsli)
            save(dfri, file=file.path(fo, paste("bsl_", i, ".RData", sep="")))
            write.table(dfri, file=file.path(fo, paste("bsl_", i, ".csv", sep="")), 
                                                    sep=";", row.names=FALSE)
            
            # save cum-merged site-table for bam[i]
            bsl <- merge(bsl, bsli)
            dfr <- as.data.frame(bsl)
            
            save(dfr, file=file.path(fo, paste("bsl_c_", i, ".RData", sep="")))
            
            write.table(dfr, file=file.path(fo, paste("bsl_c_", i, 
                    ".csv", sep="")), sep=";", row.names=FALSE)
            
            cat("\r[readSepGapTables] i:(", format(i, width=2),
                "/", n, ")\tnr sites: ",  format(size(bsl), big.mark=bm,
                                                width=9), "\n", sep="")
        }
    }
    cat("[readSepGapTables] Finished.")
}

#  Unexported and undocumented
readAccGapTables <- function(bam, profo, defo="sep_gap",
                                            idx=paste(bam, ".bai", sep=""))
{
    # setup
    require(rbamtools)
    fo <- file.path(profo, defo)
    
    if(!file.exists(fo))
        dir.create(fo)
    
    bm <- Sys.localeconv()[7]
    
    n <- length(bam)
    res <- data.frame(i=1:n, sites=numeric(n), acc=numeric(n), nov=numeric(n))
    
    for(i in 1:n)
    {
        cat("[readAccGapTables] i:(", format(i, width=2), "/", n, ")", sep="")
        
        if(!file.exists(bam[i]))
            stop("[readAccGapTables] i:", i, " File does not exist!")
        
        reader <- bamReader(bam[i])
        if(!file.exists(idx[i]))
            createIndex(reader, idx[i])
        loadIndex(reader, idx[i])
        
        if(i==1) # first bam file
        {
            bsl <- bamGapList(reader)
            dfr <- as.data.frame(bsl)
            save(dfr, file=file.path(fo, paste("bsl_", i, ".RData", sep="")))
            
            # write report values
            res$sites[1] <- dim(dfr)[1]
            res$acc[1] <- res$sites[1]
            res$nov[1] <- res$sites[1]
            
            # printout status line
            cat("\r[readAccGapTables] i:(", format(i, width=2),
                "/", n, ")\tnr sites: ",  format(size(bsl), big.mark=bm,
                                                width=9), "\n", sep="")
            
        }else{    # subsequent bam file
            
            # read sites from bam[i]
            bsli <- bamGapList(reader)
            dfri <- as.data.frame(bsli)
            
            # extract novel sites as difference from accumulated sites
            mrg <- merge(dfri[, c("id", "seqid", "lend", "rstart")], 
                        dfr[, c("id", "seqid", "lend", "rstart")], 
                        by=c("seqid", "lend", "rstart"), all.x=TRUE)
            
            mrg$new <- is.na(mrg$id.y)
            mrg$id.x <- NULL
            mrg$id.y <- NULL
            nov <- merge(dfri, mrg, all=T)
            nov <- nov[nov$new, c(4, 1, 5, 2, 3, 6:13)]
            nov$new <- NULL
            nov <- nov[order(nov$seqid, nov$lend, nov$rstart), ]
            
            # Save image
            save(dfri, nov, dfr , file=file.path(fo, paste("acc_", i,
                                                        ".RData", sep="")))
            
            # create new accumulation via merging
            bsl <- merge(bsl, bsli)
            dfr <- as.data.frame(bsl)
            
            # write report values
            res$sites[i] <- dim(dfri)[1]
            res$acc[i] <- dim(dfr)[1]
            res$nov[i] <- dim(nov)[1]
            
            # print-out status line
            cat("\r[readAccGapTables] i:(", format(i, width=2), "/", n,
                ")\tnr sites: ", format(size(bsl), big.mark=bm, width=9),
                                "\n", sep="")
        }
    }
    
    # save final image
    save(dfr, file=file.path(fo, "bsl_acc_final.RData"))
    
    cat("[readAccGapTables] Final sites: ", 
                            format(size(bsl), big.mark=bm, width=9), "\n")
    
    return(res)
}

#  Unexported and undocumented
copy_fastq <- function(infile, outfile, which, append=FALSE)
{
    if(!is.character(infile))
        stop("[copy_fastq] infile must be character!")
    
    if(!is.character(outfile))
        stop("[copy_fastq] outfile must be character!")
    
    if(!is.numeric(which))
        stop("[copy_fastq] which must be numeric!")
    
    if(any(which<=0))
        stop("[copy_fastq] Only positive numbers in which allowed!")
    
    which <- as.integer(sort(unique(which)))
    if(!is.logical(append))
        stop("[copy_fastq] append must be logical!")
    
    if(!file.exists(infile))
        stop("[copy_fastq] infile does not exist!")
    
    ans <- .Call("copy_fastq_records", infile, outfile, which,
                                            append, PACKAGE="rbamtools")
    
    bm <- Sys.localeconv()[7]
    
    if(length(which)<ans){
        cat("[copy_fastq] Incomplete copy: ", format(ans, big.mark=bm), "/", 
                format(length(which), big.mark=bm), ". EOF reached?", sep="")
    }
    return(invisible(ans))
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Function declaration for extraction of regions into separate BAM files:
#
# Extracts alignments from given (genetic) ranges and BAM files
# into a set of output BAM files.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

extractBamRegions <- function(bamFiles, ranges,
                idxFiles=paste(bamFiles, "bai", sep="."),
                outFiles=paste("bam", 1:length(bamFiles),".bam", sep=""))
{
    if(!is.character(bamFiles))
        stop("'bamFiles' must be character!")
    
    if(!is.character(idxFiles))
        stop("'idxFiles' must be character!")
    
    if(!is.character(outFiles))
        stop("'outFiles' must be character!")
    
    if(!all(file.exists(bamFiles)))
        stop("bam file(s) not found!")
    
    if(!all(file.exists(idxFiles)))
        stop("idx file(s) not found!")
    
    if(length(bamFiles) != length(idxFiles))
        stop("bamFiles and idxFiles must have equal length!")
    
    if(length(bamFiles) != length(outFiles))
        stop("bamFiles and outFiles must have equal length!")
    
    if(!all(ranges$end > ranges$start))
        stop("All 'end' positions must be greater than 'start'!")
    
    
    nr <- nrow(ranges)
    nf <- length(bamFiles)
    
    for(i in 1:nf)
    {
        reader <- bamReader(bamFiles[i], idx=TRUE)
        cat("[i] ", format(i, width=2), "\n")
        header <- getHeader(reader)
        writer <- bamWriter(header, outFiles[i])
        
        for(j in 1:nr)
        {
            sn <- getRefCoords(reader,as.character(ranges$seqid[j]))
            coords <- c(sn[1], ranges$start[j], ranges$end[j])
            brg <- bamRange(reader, coords)
            bamSave(writer, brg, refid=sn[1])        
        }
        bamClose(reader)
        bamClose(writer)
    }
    return(invisible())
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# GenomePartition class
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

.GenomePartition <- setClass("GenomePartition",
                             representation(
                                 ev="environment"),
                             prototype=prototype(
                                 ev=new.env()
                             )
)

setMethod("show", "GenomePartition", function(object)
{
    if(!exists("reftable", where=object@ev, inherits=FALSE))
        cat("An empty object of class '", class(object), "'.\n", sep="")
    else
    {
        n<-min(nrow(object@ev$reftable), 6L)
        bm<-Sys.localeconv()[7]
        cat("Object of class '", class(object), "' with ",
            format(nrow(object@ev$reftable), big.mark = bm), " rows and ",
            ncol(object@ev$reftable), " columns.\n", sep="")
        
        print(object@ev$reftable[1:n, ])
    }
})

setMethod("getRefData", "GenomePartition", function(object) {
    return(object@ev$reftable)
})

setMethod("getSeqNr", "GenomePartition", function(object) {
    return(nrow(object@ev$reftable))
})


#setGeneric("genomePartition", function(object, genome, dist=10000) standardGeneric("genomePartition"))
setMethod("genomePartition", c("bamReader","data.frame"), function(object, genome, dist=10000){
    
    colnames <- c("id", "seqid", "begin", "end")
    mtc <- match(colnames, names(genome))
    if(any(is.na(mtc)))
        stop("genome must contain columns 'id', 'seqid', 'begin', 'end'!")
    
    genome$id <- as.integer(genome$id)
    genome$begin <- as.integer(genome$begin)
    genome$end <- as.integer(genome$end)
    
    rfd <- getRefData(object)
    clv <- levels(genome$seqid)
    mtc <- match(rfd$SN, clv)
    
    if(any(is.na(mtc)))
    {
        cat("[genomePartition] Missing reader seqid in genome annotation data!\n")
        rfd <- rfd[!is.na(mtc), ]
    }
    
    reflist <- list()
    # number of chromosomes
    nc <- nrow(rfd) 
    for(i in 1:nc)
    {
        cat("[Partition] i: (", format(i, w=3), "/", nc, ") chr: ",clv[i], "\n",sep="")
        
        # Calculate Segments
        urc <- genome[genome$seqid==clv[i], ]
        coords <- as.integer(c(1, rfd$LN[i], dist, 1))
        reflist[[i]] <- .Call("get_partition_segments", urc$id,
                                urc$begin, urc$end, coords, PACKAGE="rbamtools")
    }
    
    res <- .GenomePartition()  
    rfd$nseg <- unlist(lapply(reflist, nrow))
    res@ev$reftable <- rfd
    res@ev$reflist  <- reflist
    res@ev$genome <- genome
    return(res)
})

# src = source
# setGeneric("countPartition", function(partition, src) standardGeneric("countPartition"))


# Counts on a genome partition on a single BAM file
setMethod("countPartition", c("GenomePartition", "bamReader"), function(partition, src){
    
    if(!isOpen(src))
        stop("BamReader must be opened!")
    
    # Used for formatting of console output
    bm <- Sys.localeconv()[7]
    r <- "right"
    w <- 14
    sample_name <- "sample"
    
    # Reference sequence data
    rfd <- getRefData(partition)
    nc <- getSeqNr(partition)
    
    for(i in 1:nc)
    {
        cat("[countPartition] i: (", format(i, w=3), "/", nc, ") chr: ",rfd$SN[i],sep="")
        
        # Count aligns for each partition
        ccrd <- c(refid=rfd$ID[i], begin=0, end=rfd$LN[i])
        cnt <- rangeSegCount(src, ccrd, partition@ev$reflist[[i]]$position)
        
        # Save result to partition object
        dfr<-as.data.frame(cnt)
        partition@ev$reflist[[i]][,sample_name] <- cnt@count
        partition@ev$reftable$count[i] <- sum(cnt@count)
        
        cat("\tAligns:", format(format(partition@ev$reftable$count[i], big.m=bm), w=w, j=r), "\n")
    }
    
    partition@ev$filetable <- data.frame(filename=filename(src), 
                        sample=sample_name, stringsAsFactors=FALSE)
    return(invisible())
})

setGeneric("checkPartition", function(partition, src, verbose=TRUE) standardGeneric("checkPartition"))
setMethod("checkPartition", c("GenomePartition", "data.frame"), function(partition, src, verbose=TRUE){
    
    col_names <- c("filename", "sample")
    if(any(is.na(match(col_names, names(src)))))
        stop("source data.frame must contain columns 'filename' and 'sample'!")
    
    src$filename <- as.character(src$filename)
    src$sample <- as.character(src$sample)
    
    if(any(!file.exists(src$filename)))
        stop("Missing input BAM files!")
    
    if(!all(is.na(match(src$sample,  names(partition@ev$reftable)))))
        stop("Sample names must not already be present in partition reftable!")
    
    if(!is.logical(verbose))
        stop("'verbose' must be logical")
    
    verbose <- verbose[1]
    
    rfd <- getRefData(partition)
    nf <- nrow(src)             # Number of input files
    for(i in 1:nf)
    {
        # Reference sequence data
        reader <- bamReader(src$filename[i], idx=TRUE)
        rref <- getRefData(reader)
        mtc <- match(rref$SN, rfd$SN)
        if(any(is.na(mtc)) && verbose)
        {
            message("Missing bamReader reference in partition for i=", i, ":")
            print(rref[is.na(mtc), ])
        }
        
        mtc <- match(rfd$SN, rref$SN)
        if(any(is.na(mtc)) && verbose)
        {
            message("Missing partition reference in bamReader for i=", i, ":")
            print(rfd[is.na(mtc), ])
        }
    }
})


# Counts on a Genome partition from multiple BAM files
setMethod("countPartition", c("GenomePartition", "data.frame"), function(partition, src){
    
    checkPartition(partition, src, FALSE)
    
    # Used for formatting of console output
    bm <- Sys.localeconv()[7]
    r <- "right"
    w <- 14
    
    # Reference sequence data
    rfd <- getRefData(partition)
    
    nc <- getSeqNr(partition)   # Number of Reference sequences
    nf <- nrow(src)             # Number of input files
    
    # Create empty columns in result tables
    for(i in 1:nf)
    {
        rfd[,src$sample[i]] <- 0
        partition@ev$reflist[[i]][, src$sample[i]] <- 0
    }
    
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    # Cyle all BAM files (given in src)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    for(i in 1:nf)
    {
        cat("[countPartition] file: (", 
                                format(i, width=3),"/", nf, ")\n", sep="")
        
        # Open BAM file for reading
        reader <- bamReader(src$filename[i], idx=TRUE)
        
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
        # Prepare Reference sequence (chromosome) data for opened reader:
        # Find mathes for partition sequence in BAM - Header
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
        rref <- getRefData(reader)
        mtc <- match(rfd$SN, rref$SN)
        # May contain NA's!
        rfd$rd_id <- rref$ID[mtc]
        
        
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
        # Count alignments for each sequence (chromosome) partition
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
        for(j in 1:nc)
        {    
            if(!is.na(rfd$rd_id[j]))
            {
                # cat("[countPartition] i   : (", 
                #   format(j, width=3), "/", nc, ") seq: ", 
                #   rfd$SN[j], "\r", sep="")
                
                # Coordinates for reading partition from BAM 
                ccrd <- c(refid=rfd$rd_id[j], begin=0, end=rfd$LN[j])
                
                # Do count aligns
                cnt <- rangeSegCount(reader, ccrd, 
                                        partition@ev$reflist[[j]]$position)
                dfr<-as.data.frame(cnt)
                
                # Add count values for i'th sample to partition table
                partition@ev$reflist[[j]][,src$sample[i]] <- cnt@count
                
                # Add total align count for i'th sample to reference sequence
                # table
                rfd[j, src$sample[i]] <- sum(cnt@count)
            }
        }
    }
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    # Add data to incoming partition object (and return nothing)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    partition@ev$reftable <- rfd
    partition@ev$filetable <- src
    cat("[countPartition] Finished.\n")
    return(invisible())
})


#setGeneric("getFileTable", function(object) standardGeneric("getFileTable"))
setMethod("getFileTable", "GenomePartition", function(object) 
    { return(object@ev$filetable) })


# setGeneric("getAlignCounts", function(object) standardGeneric("getAlignCounts"))
setMethod("getAlignCounts","GenomePartition", function(object){
    
    if(is.null(object@ev$filetable))
        stop("No BAM aligns counted yet!")
    
    # Number of bam files
    filetable <- object@ev$filetable
    nf <- nrow(filetable)
    
    nc <- getSeqNr(object)   # Number of Reference sequences
    
    # Create output columns in genome table
    for(j in 1:nf)
        object@ev$genome[, filetable$sample[j]] <- 0
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    # !! id values in returned table correspond to  !!
    # !! source_id values in reflist                !!
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    for(i in 1:nc)      # For each chromosome
    {
        cnt <- object@ev$reflist[[i]][object@ev$reflist[[i]]$source != 0, ]
        mtc <- match(cnt$source_id, object@ev$genome$id)
        
        for(j in 1:nf)  # Copy one count column for each file
            object@ev$genome[, filetable$sample[j]][mtc] <- cnt[, filetable$sample[j]]
    }
    return(invisible(object@ev$genome))
})


setMethod("getGridAlignCounts", "GenomePartition", function(object){
    nc <- getSeqNr(object)
    l <- list()
    for(i in 1:nc)
        l[[i]] <- object@ev$reflist[[i]][object@ev$reflist[[i]]$source==0, ]
    
    res <- do.call(rbind, l)
    res$source <- NULL
    res$source_id <- NULL
    
    return(invisible(res))
})


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#  End of File (rbamtools.r)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
