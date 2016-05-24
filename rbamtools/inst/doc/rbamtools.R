### R code from vignette source 'rbamtools.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: rbamtools.Rnw:64-67
###################################################
library(rbamtools)
library(xtable)
options(width=60)


###################################################
### code chunk number 2: rbamtools.Rnw:591-595
###################################################
bam <- system.file("extdata", 
                "accepted_hits.bam", package="rbamtools")
# Open bam file
reader <- bamReader(bam)


###################################################
### code chunk number 3: rbamtools.Rnw:606-608 (eval = FALSE)
###################################################
## bamSort(reader, prefix="my_sorted", 
##             byName=FALSE, maxmem=1e+9)


###################################################
### code chunk number 4: rbamtools.Rnw:615-616 (eval = FALSE)
###################################################
## createIndex(reader, idx_filename="index_file_name.bai")


###################################################
### code chunk number 5: rbamtools.Rnw:623-624 (eval = FALSE)
###################################################
## createIndex(reader)


###################################################
### code chunk number 6: rbamtools.Rnw:631-633
###################################################
idx <- system.file("extdata", "accepted_hits.bam.bai", package="rbamtools")
loadIndex(reader, idx)


###################################################
### code chunk number 7: rbamtools.Rnw:636-637
###################################################
indexInitialized(reader)


###################################################
### code chunk number 8: rbamtools.Rnw:641-642
###################################################
reader <- bamReader(bam, idx=TRUE)


###################################################
### code chunk number 9: rbamtools.Rnw:650-651
###################################################
getRefData(reader)


###################################################
### code chunk number 10: rbamtools.Rnw:664-668 (eval = FALSE)
###################################################
## header <- getHeader(reader)
## writer <- bamWriter(header,"test.bam")
## # Write alignments using bamSave
## bamClose(writer)


###################################################
### code chunk number 11: rbamtools.Rnw:695-697
###################################################
header <- getHeader(reader)
htxt <- getHeaderText(header)


###################################################
### code chunk number 12: rbamtools.Rnw:720-739
###################################################
bh <- new("bamHeaderText")

headl <- new("headerLine")
setVal(headl, "SO", "coordinate")

dict <- new("refSeqDict")
addSeq(dict, SN="chr1",  LN=249250621)
addSeq(dict, SN="chr16", LN=90354753)
dict

prog <- new("headerProgram")
setVal(prog, "ID", "TopHat")
setVal(prog, "PN", "tophat")
setVal(prog, "CL",
    "tophat --library-type fr-unstranded hs_ucsc_index reads.fastq")
setVal(prog, "DS", "Description")
setVal(prog, "VN", "2.0.0")
bh <- bamHeaderText(head=headl, dict=dict, prog=prog)
header <- bamHeader(bh)


###################################################
### code chunk number 13: rbamtools.Rnw:747-748
###################################################
align <- getNextAlign(reader)


###################################################
### code chunk number 14: rbamtools.Rnw:777-788 (eval = FALSE)
###################################################
## name(align)
## flag(align)
## refID(align)
## position(align)
## mapQuality(align)
## cigarData(align)
## nCigar(align)
## mateRefID(align)
## matePosition(align)
## alignSeq(align)
## alignQual(align)


###################################################
### code chunk number 15: rbamtools.Rnw:815-826 (eval = FALSE)
###################################################
## paired(align)
## properPair(align)
## unmapped(align)
## mateUnmapped(align)
## reverseStrand(align)
## mateReverseStrand(align)
## firstInPair(align)
## secondInPair(align)
## secondaryAlign(align)
## failedQC(align)
## pcrORopt_duplicate(align)


###################################################
### code chunk number 16: rbamtools.Rnw:831-832
###################################################
unmapped(align) <- TRUE


###################################################
### code chunk number 17: rbamtools.Rnw:841-850
###################################################
align <- bamAlign("HWUSI-0001", "ATGTACGTCG", "Qual/Strng",
                "4M10N6M", refid=0, position=100)
align
name(align)
alignSeq(align)
alignQual(align)
cigarData(align)
refID(align)
position(align)


###################################################
### code chunk number 18: rbamtools.Rnw:886-890
###################################################
coords <- c(0, 899000, 900000)
names(coords) <- c("refid","start","stop")
range <- bamRange(reader,coords)
size(range)


###################################################
### code chunk number 19: rbamtools.Rnw:897-902
###################################################
getRefData(reader)
coords <- c(0,0,249250621)
names(coords) <- c("refid","start","stop")
range <- bamRange(reader,coords)
size(range)


###################################################
### code chunk number 20: rbamtools.Rnw:907-911
###################################################
coords <- getRefCoords(reader,"chr1")
coords
range <- bamRange(reader,coords)
size(range)


###################################################
### code chunk number 21: rbamtools.Rnw:921-926
###################################################
range
getCoords(range)
getSeqLen(range)
getParams(range)
getRefName(range)


###################################################
### code chunk number 22: rbamtools.Rnw:932-933
###################################################
getAlignRange(range)


###################################################
### code chunk number 23: rbamtools.Rnw:952-953
###################################################
align <- getNextAlign(range)


###################################################
### code chunk number 24: rbamtools.Rnw:960-966 (eval = FALSE)
###################################################
## rewind(range)
## while(!is.null(align))
## {
##   # Process align data here
##   align <- getNextAlign(range)
## }


###################################################
### code chunk number 25: rbamtools.Rnw:972-973
###################################################
rdf <- as.data.frame(range)


###################################################
### code chunk number 26: rbamtools.Rnw:994-999
###################################################
coords <- getRefCoords(reader, "chr1")
gl <- gapList(reader, coords)
gl
dfr <- as.data.frame(gl)
dfr[1:6, c(1:3, 5:8)]


###################################################
### code chunk number 27: rbamtools.Rnw:1009-1012 (eval = FALSE)
###################################################
## size(gl)
## nAligns(gl)
## nAlignGaps(gl)


###################################################
### code chunk number 28: rbamtools.Rnw:1055-1063
###################################################
coords <- getRefCoords(reader, "chr1")
sl <- siteList(reader, coords)
size(sl)
nAligns(sl)
nAlignGaps(sl)
sl
df <- as.data.frame(sl)
head(df)


###################################################
### code chunk number 29: rbamtools.Rnw:1077-1085
###################################################
bsl <- bamGapList(reader)
bsl
size(bsl)
nAligns(bsl)
nAlignGaps(bsl)
summary(bsl)
dfr <- as.data.frame(bsl)
head(dfr)


###################################################
### code chunk number 30: rbamtools.Rnw:1099-1102
###################################################
bam<-system.file("extdata","accepted_hits.bam",package="rbamtools")
rpb<-readPooledBamGaps(bam)
rpdf<-readPooledBamGapDf(bam)


###################################################
### code chunk number 31: rbamtools.Rnw:1112-1113
###################################################
xtable(head(rpdf)[, c(1:6, 8, 9, 13)])


###################################################
### code chunk number 32: hist_gqs
###################################################
plot(table(rpdf$gqs), type="h", las=1, bty="n", lwd=2, 
        xlab="gqs", ylab="Number of gap sites",
        main="Distribution of gqs values")


###################################################
### code chunk number 33: rbamtools.Rnw:1147-1233 (eval = FALSE)
###################################################
## scanGapSites <- function(bam, yieldSize=1e6, mc.cores=2)
## {
##     require(Rsamtools)
##     require(GenomicAlignments)
##     
##     mc.cores <- as.integer(mc.cores)
##     
##     # Function will be called by mclapply
##     doScanBam <- function(bamFile)
##     {
##         open(bamFile)
##         
##         # Create empty container.
##         gPos <- GRanges()
##         
##         # Fill container by processing 'yieldSize' reads at a time
##         while(
##             sum(
##                 elementLengths(
##                     records <- scanBam(
##                                 bamFile,
##                                 param=ScanBamParam(
##                                     flag=scanBamFlag(isUnmappedQuery=FALSE),
##                                     what=scanBamWhat()[c(3, 5, 8)]
##                                     )
##                                 )[[1]]
##                     )
##             ) > 0
##         ){
##             
##             nOps <- cigarRangesAlongReferenceSpace(records$cigar,ops="N")
##             sel <- elementLengths(nOps) != 0
##             gPos <- c(gPos, 
##                         GRanges(seqnames=rep(records$rname[sel], 
##                                             elementLengths(nOps)[sel]),
##                                             ranges=unlist(shift(nOps[sel], 
##                                             records$pos[sel]))
##                             )
##                     )
##         }
##         close(bamFile)
##         # Return all gap positions
##         return(gPos)
##     }
##     
##     cat("[scanGapSites] Processing", length(bam), "Files.\n")
##     
##     
##     bamFileList <- BamFileList(bam, yieldSize = yieldSize)
##     gList <- mclapply(bamFileList, doScanBam)
##     
##     sz <- object.size(gList)
##     bm<-Sys.localeconv()[7]
##     cat("[scanGapSites] Collected object of size", 
##             format(as.numeric(object.size(gList)), big.mark=bm),
##             "bytes.\n")
##     
##     # - - - - - - - - - - - - - - - - - - - - - - - - - #
##     # Get all unique positions across all samples
##     # - - - - - - - - - - - - - - - - - - - - - - - - - #
##     uPos <- unique(Reduce("c", gList))
##     
##     # - - - - - - - - - - - - - - - - - - - - - - - - - #
##     # Create the count table by 
##     # transforming the ranges into character strings.
##     # - - - - - - - - - - - - - - - - - - - - - - - - - #
##     ref <- paste(seqnames(uPos), start(uPos), end(uPos), sep="-")
##     
##     # Will be called by mclapply
##     doTable <- function(grng, ref)
##     {
##         tab <- table(paste(seqnames(grng), start(grng), end(grng),sep="-"))
##         tab[match(ref,names(tab))]
##     }
##     
##     count.table <- do.call("cbind", 
##                     mclapply(gList, doTable, ref, mc.cores=mc.cores))
##     rownames(count.table) <- ref
##     
##     cat("[scanGapSites] Number of sites:", 
##             format(nrow(count.table), big.mark=bm),
##             ".\n")
##     
##     cat("[scanGapSites] Finished.\n")
##     return(count.table)
## }


###################################################
### code chunk number 34: rbamtools.Rnw:1301-1304
###################################################
coords <- c(0, 0, 14730)
count <- bamCount(reader, coords)
xtable(matrix(count, nrow=1))


###################################################
### code chunk number 35: rbamtools.Rnw:1308-1309
###################################################
count <- bamCountAll(reader, verbose=TRUE)


###################################################
### code chunk number 36: rbamtools.Rnw:1312-1313
###################################################
xtable(count, digits=0)


###################################################
### code chunk number 37: rbamtools.Rnw:1324-1327
###################################################
align <- bamAlign("HWUSI-0001", "ACCGGGTTTT","Qual/Strng",
                            "4M10N6M", refid=0, position=100)
countNucs(align)


###################################################
### code chunk number 38: rbamtools.Rnw:1330-1334
###################################################
reader <- bamReader(bam, idx=TRUE)
coords <- c(0, 0, 14730)
range <- bamRange(reader, coords)
countNucs(range)


###################################################
### code chunk number 39: rbamtools.Rnw:1349-1350
###################################################
ncs <- nucStats(reader)


###################################################
### code chunk number 40: rbamtools.Rnw:1353-1354
###################################################
xtable(ncs, digits=c(0, 0, 0, 0, 0, 0, 0, 2, 2))


###################################################
### code chunk number 41: rbamtools.Rnw:1362-1363
###################################################
ncs <- nucStats(bam)


###################################################
### code chunk number 42: rbamtools.Rnw:1366-1367
###################################################
xtable(ncs, digits=c(0, 0, 0, 0, 0, 0, 0, 2, 2))


###################################################
### code chunk number 43: rbamtools.Rnw:1387-1388 (eval = FALSE)
###################################################
## createIdxBatch(bam)


###################################################
### code chunk number 44: rbamtools.Rnw:1412-1419 (eval = FALSE)
###################################################
## reader <- bamReader(bam)
## readerToFastq(reader, "out.fastq")
## bamClose(reader)
## # Reopen in order to point to first alignment
## reader <- bamReader(bam)
## index <- sample(1:100, 20)
## readerToFastq(reader, "out_subset.fastq", which=index)


###################################################
### code chunk number 45: rbamtools.Rnw:1429-1435 (eval = FALSE)
###################################################
## reader <- bamReader(bam, idx=TRUE)
## coords <- as.integer(c(0,0,249250621))
## range <- bamRange(reader,coords)
## rangeToFastq(range,"rg.fq.gz")
## index <- sample(1:size(range),100)
## rangeToFastq(range,"rg_subset.fq.gz",which=index)


###################################################
### code chunk number 46: rbamtools.Rnw:1451-1456
###################################################
qdf <- getQualDf(range)
qdf[32:38,1:10]
qdr <- getQualDf(range,prob=TRUE)
qrr <- round(qdr,2)
qrr[32:38,1:10]


###################################################
### code chunk number 47: rbamtools.Rnw:1463-1465
###################################################
qt <- getQualQuantiles(range,c(0.25,0.5,0.75))
qt[,1:10]


###################################################
### code chunk number 48: rbamtools.Rnw:1471-1472
###################################################
plotQualQuant(range)


###################################################
### code chunk number 49: rbamtools.Rnw:1492-1509
###################################################
# WASH7P coordinates
xlim <- c(10000, 30000)
coords <- c(0,xlim[1], xlim[2])
range <- bamRange(reader, coords)
bamClose(reader)
ad <- alignDepth(range)
ad
getParams(ad)
# Identifier
gene <- "WASH7P"
ensg_id <- "ENSG00000227232"
enst_id <- "ENST00000538476"
# Get exon positions
start <- c(14411, 15000, 15796, 15904, 16607, 16748, 16858, 17233,
                    17602, 17915, 18268, 24737, 29534)
end <-   c(14502, 15038, 15901, 15947, 16745, 16765, 17055, 17364,
                    17742, 18061, 18366, 24891, 29806)


###################################################
### code chunk number 50: rbamtools.Rnw:1512-1520
###################################################
plotAlignDepth(ad, lwd = 2, xlim = xlim,
            main = paste("Align depth for gene",gene),
            ylab = "Align depth", start = start,
            end = end, strand = "-",
            transcript = paste("Chromosome 1",
                "\tGene ENSG00000227232", ensg_id, 
                "\tTranscript ",enst_id
))


###################################################
### code chunk number 51: rbamtools.Rnw:1550-1561
###################################################
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# B) Count range segment
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
reader <- bamReader(bam, idx=TRUE)
coords <- c(0, 0, 2e4)
segments <- seq(14000, 20000, 20)
segcount<-rangeSegCount(reader, coords, segments)
segcount
dfr<-as.data.frame(segcount)
sum(dfr$count)



###################################################
### code chunk number 52: rbamtools.Rnw:1564-1569
###################################################
plot(count~position, dfr, type="l", 
        las=1, bty="n", lwd=1.5, col="dodgerblue2",
        xlab="Position on Chromosome 1",
        ylab="Alignment count",
        main="Number of alignments in genomic segments of 20 nucleotides size")


