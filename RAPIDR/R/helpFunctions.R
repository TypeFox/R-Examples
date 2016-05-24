#require(Rsamtools)
#require(GenomicRanges)
#require(BSgenome.Hsapiens.UCSC.hg19)

#' @title makeGCContentData 
#' 
#' @description
#' This function calculates the GC content using the hg19 
#' reference genome in bins of user-defined size 
#' 
#' @param gc.fname file name of the output file. The file will 
#'        be written as a .RData object 
#' @param k bin size, default is 20,000 bp 
#' @export 
#' 
makeGCContentData <- function (gc.fname, k = 20000) { 
  
  library(BSgenome.Hsapiens.UCSC.hg19)
  hg19 <- BSgenome.Hsapiens.UCSC.hg19 
  chrnames <- seqnames(hg19)[1:24]
  gcContent <- list() 
  for (c in chrnames) { 
    message("Processing for ", c)
    thisgcContent <- letterFrequencyInSlidingView(getSeq(hg19, c), k, "GC", as.prob=TRUE)
    chrLen <- seqlengths(hg19)[c]
    pos <- seq(1,chrLen,k)
    gcContent[[c]] <- thisgcContent[pos]
  }
 
  save(gcContent, file=gc.fname)
  return(gcContent)
}

#' @title makeGCContentPerPos 
#' 
#' @description
#' This function calculates the GC content using the hg19 
#' reference genome for each position in the chromosome  
#' 
#' @param gc.fname file name of the output file. The file will 
#'        be written as a .RData object 
#' @param k bin size, default is 200 bp 
#' @importFrom Biostrings letterFrequencyInSlidingView
#' @import GenomicRanges 
#' 
#' 
makeGCContentPerPos <- function (gc.fname, k = 200) { 
  
  library(BSgenome.Hsapiens.UCSC.hg19)
  hg19 <- BSgenome.Hsapiens.UCSC.hg19 
  chrnames <- seqnames(hg19)[1:24]
  #gcContent <- list() 
  for (c in chrnames) { 
     message("Processing for ", c)
     thisgcContent <- letterFrequencyInSlidingView(getSeq(hg19, c), k, "GC") 
     fname <- paste(gc.fname, c, ".Rdata", sep = "") 
     save(thisgcContent, file = fname) 
     #gcContent[[c]] <- thisgcContent
  }
 
}


#############################################################
## Internal functions - not exported 
#############################################################

find.chrY.excl.bins <- function (cleaned.good.ratios, bin.names, cleaned.sampleIDs.with.outcomes ) {
  chrY.pos <- which(bin.names == "chrY")
  not.chrY.pos <- which(bin.names != "chrY")
  females  <- which(cleaned.sampleIDs.with.outcomes$Gender == "Female" & 
                     cleaned.sampleIDs.with.outcomes$Dx == "Normal") 
 
  female.cts.mean <- colMeans(cleaned.good.ratios[females, ])
  female.cts.mean[not.chrY.pos] <- 0 
  mean.cts.asc <- order(-female.cts.mean)
  per_cleaned <- sum(female.cts.mean[mean.cts.asc[1:25]])*100/sum(female.cts.mean[mean.cts.asc])
  message("... excluding 25 bins,  ", round(per_cleaned, digits=1), " % of counts mapped to chrY in females." )
  return(mean.cts.asc[1:25])
}

find.excl.bins <- function (cleaned.good.ratios, bin.names, cleaned.sampleIDs.with.outcomes ) {
  not.chrXY.pos <- which(!bin.names %in% c("chrY", "chrX"))
  cts.mean <- colMeans(cleaned.good.ratios)  
  cts.sd   <- apply(cleaned.good.ratios, 2, sd) 

  overall.cts.sd <- sd(cts.mean) 
  cts.median <- median(cts.mean) 
  #message("median count is: ", cts.median) 
  #message("sd count is: ", overall.cts.sd) 
  excess.pos <- which(cts.mean > cts.median + overall.cts.sd) 
  cts.mean[excess.pos] <- 0 
  cts.sd[excess.pos] <- 0 
 
  real.sd <- sqrt(cts.mean*length(cts.mean)/1e6) 

  sd.of.sd <- sd(real.sd) 
  #message("sd limit: ", sd.of.sd) 
  excess.sd.pos <- which(cts.sd > real.sd + sd.of.sd) 
  #message(length(excess.sd.pos), " over variance")  
  cts.mean[excess.sd.pos] <- 0 
  cts.sd[excess.sd.pos] <- 0 

  #pdf("plot_mean_sd_oldgc.pdf") 
  #plot(cts.mean, cts.sd*cts.sd) 
  #dev.off()  

  excl.bins <- unique(c(excess.sd.pos, excess.pos) )  
  excl.bins <- subset(excl.bins, excl.bins %in% not.chrXY.pos) 
  message(length(excl.bins), " bins excluded altogether")  
  chr13.excl.bins <- which(bin.names[excl.bins] == "chr13") 
  chr21.excl.bins <- which(bin.names[excl.bins] == "chr21") 
  chr18.excl.bins <- which(bin.names[excl.bins] == "chr18") 
  message("Bins excluded from chr13: ", length(chr13.excl.bins)) 
  message("Bins excluded from chr18: ", length(chr18.excl.bins)) 
  message("Bins excluded from chr21: ", length(chr21.excl.bins)) 
  
  return(excl.bins)   
}
 
#' @title ReadBed 
#' 
#' @description
#' This function loads a bed file 
#' 
#' @param bedfile name of the bed file 
#' @import GenomicRanges
#' 
ReadBed <- function ( bedfile ) {
  data <- read.table(bedfile, header=FALSE) 
  colnames(data) <- c('chr','start','end') 
  bed <- with(data, GRanges(chr, IRanges(start, end))) 
  return(bed)
}

#' @title BinBam
#' 
#' @description
#' Given a list of bam files, this function writes the output 
#' to a text file after binning and doing gc Correction 
#' 
#' @param bam.file name of bamfile 
#' @param index index of bam file 
#' @param mask  mask file in the bed format 
#' @param k bin size in kilobases 
#' @import GenomicRanges Rsamtools GenomicAlignments
#' 
BinBam <- function(bam.file, index = bam.file, mask = NULL, k = 20000) {  

  library(BSgenome.Hsapiens.UCSC.hg19)
  hg19 <- BSgenome.Hsapiens.UCSC.hg19 
  bf   <- BamFileList(bam.file, index = index)
   
  chr.lens <- seqlengths(hg19)[1:24]
  tiles<- tileGenome(chr.lens, tilewidth = k, cut.last.tile.in.chrom=TRUE)

  message("doing the binning")
  olap <- summarizeOverlaps(tiles, bf)
  counts.mat <- assays(olap)$counts
  chrnames <- as.character(seqnames(tiles))
  countData <- data.frame(counts = counts.mat[,1], mask.counts = counts.mat[,1],chrname = chrnames ) 
  
  if ( !is.null(mask) ) {
    
     message("Binning using the mask file")
     bed.granges <- ReadBed(mask)
     all.granges <- readGAlignments(bam.file)    
     overlaps    <- subsetByOverlaps(all.granges, bed.granges)
     olap        <- countOverlaps(tiles, overlaps)
     countData$mask.counts <-  countData$counts - olap
     
  } 
  return(countData)
   
}

BinListOfBam <- function (bam.file.list, mask = NULL, index = bam.file.list, k = 20000) {
  
  all.binnedCounts <- NULL 

  for (each.bam in bam.file.list) { 
     ptm <- proc.time()
     countData <- BinBam(each.bam, mask = mask, k = k) 
     gc() 
     if (is.null(all.binnedCounts)) { 
       all.binnedCounts <- data.frame(t(countData$counts))
       masked.binnedCounts <- data.frame(t(countData$mask.counts))
     } else { 
       all.binnedCounts <- rbind(all.binnedCounts, countData$counts) 
       masked.binnedCounts <- rbind(masked.binnedCounts, countData$mask.counts)
     }

     names(all.binnedCounts) <- countData$chrname
     names(masked.binnedCounts) <- countData$chrname
     elapsed <- (proc.time() - ptm)["elapsed"]
     message("Binning done in ", elapsed)    
  }
  
  return(list(all.binnedCounts, masked.binnedCounts))
}

aggregateCounts <- function (all.binnedCounts) {
  bam.names <- subset(names(all.binnedCounts),names(all.binnedCounts) != "chrname")
  aggregate.counts <- data.frame(chrname = levels(all.binnedCounts$chrname))
  
  for (each.bam in bam.names) {
    tmp.results <- aggregate(all.binnedCounts[each.bam], list(chrname = all.binnedCounts$chrname), sum)
    message("Next...", each.bam)
    aggregate.counts <- merge(tmp.results, aggregate.counts, by = "chrname")
  }
  
  return(aggregate.counts)
  
}

gcCorrectCounts <- function (count.data, bin.names, masked.count.data, gcContentFile) { 
  
  load(gcContentFile)

  if (length(gcContent) != 24) { 
    stop("gcContent file not of the right format. Try running makeGCContentData()")
  }
  
  if (sum(sapply(gcContent,length)) != ncol(count.data)) {
    stop("Number of bins in gcContent file does not equal the number of bins in count.data")    
  }
  
  gcContent <- unlist(gcContent)
  
  output.counts <- count.data
  nsamples <- nrow(count.data) 

  for ( i in 1:nsamples ) {
    message("GC correcting sample: ",i)
    counts.with.gc <- data.frame(chrname = bin.names, counts = count.data[i,], masked.counts = masked.count.data[i,], gc = gcContent)
    names(counts.with.gc)[1] <- "chrname"
    output.counts[i,] <- gcCorrectOneSample(counts.with.gc)
    #output.counts[i,] <- gcCorrectOneSampleLoess(counts.with.gc)
  }

  #save(output.counts, count.data, file = "gc.Rdata")
  #load("gc.Rdata")
  return(output.counts)
}

gcCorrectOneSample <- function (counts.with.gc) {
  chr <- rep(0, 22)
  for (i in c(1:22)){
    chr[i] <-paste("chr", i, sep="") 
  } 
  
  brks<-seq(0,1,0.005) 
  #print(counts.with.gc)
  counts.with.gc$bin <- findInterval(counts.with.gc$gc, brks)
  
  # Filter the bins with GC content of 0 and 
  # anything other than chr1-22  

  sub.counts.with.gc<-subset(counts.with.gc,counts.with.gc$chrname %in% chr & counts.with.gc$gc > 0)  
  # Calculate the mean of the counts by aggregating by bin 
  # Filter out the outliers in each bin 
  counts.avg <- mean(sub.counts.with.gc$counts)
  counts.sd  <- sd(sub.counts.with.gc$counts)

  t <- tempfile()
  pdf(file = t)
  outliers <- which(sub.counts.with.gc$counts %in% boxplot(sub.counts.with.gc$counts)$out)
  dev.off() 

  sub.counts.with.gc <- sub.counts.with.gc[-outliers,]
  gcCorrect.tab <- aggregate(counts ~ bin, data = sub.counts.with.gc, FUN = mean, na.action = NULL )
  gcCorrect.tab$gcWeights<- counts.avg/gcCorrect.tab$counts
  # Clean up the Inf values in case there are no counts in that GC bin 
  gcCorrect.tab$gcWeights[! is.finite(gcCorrect.tab$gcWeights)] <- 1.0 

  counts.with.gc$id <- 1:nrow(counts.with.gc)
  new.counts<-merge(counts.with.gc, gcCorrect.tab, all.x = TRUE, by="bin")
  #new.counts$corrected.counts <- new.counts$counts.x * new.counts$gcWeights 
  new.counts$corrected.counts <- new.counts$masked.counts * new.counts$gcWeights
  new.counts$corrected.counts[is.na(new.counts$corrected.counts)] <- 0.0 
  new.counts$gcWeights[is.na(new.counts$gcWeights)] <- 0.0 
  new.counts <- new.counts[order(new.counts$id),] 
  return(new.counts$corrected.counts)
}

gcCorrectOneSampleLoess <- function (counts.with.gc) {
  chr <- rep(0, 22)
  for (i in c(1:22)){
    chr[i] <-paste("chr", i, sep="") 
  } 
  
  print(counts.with.gc$counts[100:150])
  brks<-seq(0,1,0.005) 
  #print(counts.with.gc)
  counts.with.gc$bin <- findInterval(counts.with.gc$gc, brks)

  
  # Filter the bins with GC content of 0 and 
  # anything other than chr1-22  
  
  sub.counts.with.gc<-subset(counts.with.gc,counts.with.gc$chrname %in% chr & counts.with.gc$gc > 0)  
  # Calculate the mean of the counts by aggregating by bin 
  # Filter out the outliers in each bin 
  counts.avg <- mean(sub.counts.with.gc$counts)
  outliers <- which(sub.counts.with.gc$counts %in% boxplot(sub.counts.with.gc$counts)$out)
  
  sub.counts.with.gc <- sub.counts.with.gc[-outliers,]
  ptm <- proc.time() 
  #lw <- loess(counts ~ bin, data = sub.counts.with.gc)
  #regress.counts <- data.frame(regressed.counts = predict(lw, 1:200), bin = 1:200)
  lw <- lowess(sub.counts.with.gc$bin, sub.counts.with.gc$counts, delta = 0.1)
  regress.df <- data.frame(bin = lw[[1]], regressed.counts = lw[[2]])
  regress.counts <- regress.df[!duplicated(regress.df),] 
  print(proc.time() - ptm)

  regress.counts$gcWeights<- counts.avg/regress.counts$regressed.counts
  counts.with.gc$id <- 1:nrow(counts.with.gc)
  new.counts <- merge(counts.with.gc, regress.counts, all.x = TRUE, by = "bin")
  new.counts$corrected.counts <- new.counts$masked.counts * new.counts$gcWeights
  #new.counts$corrected.counts <- counts.avg + (new.counts$masked.counts - new.counts$regressed.counts) 
  new.counts$corrected.counts[is.na(new.counts$corrected.counts)] <- 0.0 
  new.counts <- new.counts[order(new.counts$id),] 
  print(new.counts$corrected.counts[100:150])
  return(new.counts$corrected.counts)
}
  
