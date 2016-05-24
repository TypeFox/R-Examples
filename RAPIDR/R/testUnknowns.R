#' @title testUnknowns
#' 
#' @description
#' This function uses a reference set to test the unknown samples 
#' The unknown samples needs to be a binned counts file which can 
#' be created using the \link{makeBinnedCountsFile} function.  
#' 
#' @param ref.data.set rapidr.ref object which contains the baselines and 
#'        the corrections used to create the baselines  
#' @param unknowns.counts.file file name of the file with the binned counts 
#'        of the unknowns, first column needs to be the sampleID, first row 
#'        needs to be the chromosome names of each bin  
#' @param gcContentFile file name of a .Rdata object which contains 
#'        information on GC content in the genome       
#' @param combined.counts.fname file name to write to for the combined counts 
#'        per chromosome 
#' @param masked.counts.file optional file of the binned counts after masking 
#' 
#' @return data.frame with z-scores for chr21, chr18, chr13, and the fetal sex 
#'         which can be male, female or monosomy X. For males, there is also
#'         an estimate of the fetal fraction using the deficit of chrX  
#' 
#' @importFrom data.table fread 
#' @export 
#' 
testUnknowns <- function( ref.data.set, unknowns.counts.file, 
                          gcContentFile = NULL, 
                         masked.counts.file = NULL, combined.counts.fname = NULL ) {

  # Check that the input reference set is of the right class 
  if ( class(ref.data.set) != "rapidr.ref") { 
     stop("ref.data.set is not of the rapidr class? Check.")
  }

  PCA_output <-  ref.data.set[["PCA"]]
  baselines  <-  ref.data.set[["baselines"]]
  excl.bins  <-  ref.data.set[["excl.bins"]] 
  gcCorrect  <-  ref.data.set[["do.gcCorrect"]]
  PCA        <-  ref.data.set[["do.PCA"]]
  
  # Make sure there is a gcContentFile provided if doing gc correction 
  if (gcCorrect == TRUE) {
    if (is.null(gcContentFile)) {
      stop("You need to provide a gcContent file. Try running makeGCContentData()")
    }
  }
  
  binned.counts <- fread(unknowns.counts.file, sep=',', colClasses=list(character=1), header=TRUE)
  #n.cols <- ncol(binned.counts)
  #ind<-rep(FALSE, n.cols)
  #ind[1] <-TRUE

  bin.names <- names(binned.counts) 
  bin.names <- bin.names[2:length(bin.names)] 

  binned.counts <- data.frame(binned.counts) 
  sampleIDs <- binned.counts[,1] 
  binned.counts <- binned.counts[,-c(1)] 
 
  n.all.samples <- nrow(binned.counts)
  n.bins <- ncol(binned.counts)
  binned.counts <- as.matrix(binned.counts) 
 
  if (!is.null(masked.counts.file)) {
    masked.binned.counts <- fread(masked.counts.file, sep=',', colClasses=list(character=1), header = TRUE )
    masked.binned.counts <- data.frame(masked.binned.counts) 
    masked.binned.counts <- masked.binned.counts[,-c(1)] 
    masked.binned.counts <- sapply(masked.binned.counts, as.numeric) 
    masked.binned.counts <- as.matrix(masked.binned.counts)  
  }
  
  if ( gcCorrect & !is.null(masked.counts.file)) {
    message("Doing gc correction")
    binned.counts <- gcCorrectCounts(binned.counts, bin.names, masked.binned.counts,  gcContentFile)
  } else if ( gcCorrect & is.null(masked.counts.file)) {
    binned.counts <- gcCorrectCounts(binned.counts, bin.names, binned.counts,  gcContentFile)
  } else if (!is.null(masked.counts.file)) {
    # If there is a masked binned counts file, use that 
    binned.counts <- masked.binned.counts 
  }

  # Setting bins in the exclusion list to have zero counts 
  binned.counts[,excl.bins] <- 0.0 
  
  total.counts <-rowSums(binned.counts) 
  binned.ratios <- as.matrix(binned.counts)
  discard.pos <- c()
  n.discard <- 0  
  
  rm(binned.counts) # trying to save memory 
  
  for (i in 1:n.all.samples) { 
    this.total<-total.counts[i] 
    # Multiply by 1million to avoid small ratio numbers 
    binned.ratios[i,] <- binned.ratios[i,]/this.total * 1e6     
  }
  
  if ( PCA ) {
    cleaned.ratios <- correct.counts.with.PCA(PCA_output, binned.ratios)
    rm(binned.ratios)
  } else { 
    cleaned.ratios <- binned.ratios 
    rm(binned.ratios)
  }
  
  cleaned.counts.per.chr <- sum.counts(cleaned.ratios, bin.names, 
                                       total.counts, sampleIDs)

  message("Doing QC")
  qc.results <- data.qc(cleaned.counts.per.chr, sampleIDs, baselines)
  
  unknowns.results <- callUnknowns( cleaned.counts.per.chr, sampleIDs, baselines )

  # If there is a file name provided, write the cleaned, combined counts to a file 
  if ( ! is.null(combined.counts.fname) ) {
    write.csv(cleaned.counts.per.chr, file = combined.counts.fname, quote = FALSE, row.names = FALSE)     
  }
    
  rapidr.test <- list() 
  rapidr.test[["results"]] <- unknowns.results 
  rapidr.test[["qc"]] <- qc.results 
  rapidr.test[["baselines"]] <- baselines 
 
  class(rapidr.test) <- "rapidr.test"
  return(rapidr.test)
}

#' @title writeResultsToFile
#' 
#' @description
#' This function writes the rapidr.test object to a file 
#' 
#' @param rapidr.test rapidr.test object which contains the results from 
#'         testing a sample 
#' @param output.fname file name of the output file (it will be tab delimited)
#' 
#' @export 
#' 
writeResultsToFile <- function(rapidr.test, output.fname) {
  
  if (class(rapidr.test) != "rapidr.test") {
     stop("Input needs to be a rapidr.test object, the output from the function testUnknowns()")
  } 

  unknowns.results <- rapidr.test[["results"]]
  qc.results       <- rapidr.test[["qc"]] 

  results.tab <- cbind(unknowns.results, qc.results)
  # Write output to a file if a filename is given 
  message("Writing output to file... ")
  write.table(results.tab, file = output.fname, quote = FALSE, row.names = FALSE, sep = "\t" )

}

#' @title plotTestSample
#' 
#' @description
#' This function plots the QC information. Note that it can only plot one sample at a time 
#' 
#' @param rapidr.test rapidr.test object which contains the results from 
#'         testing a sample 
#' @param input.sampleID the sampleID that you like to plot 
#' @param ordering normal or gc, default is normal, when set to GC, the chromosomes are 
#'        ordered according to GC content
#' @export 
#' 
plotTestSample <- function(rapidr.test, input.sampleID, ordering = "normal") {

  input.sampleID <- as.character(input.sampleID)
  test.results <- rapidr.test[["results"]]
  qc.results <- rapidr.test[["qc"]]
  baselines  <- rapidr.test[["baselines"]]
  normals.ratios.mean <- baselines[["normals.mean"]]
  normals.ratios.sd   <- baselines[["normals.sd"]] 

  if (!input.sampleID %in% qc.results$sampleID) {
    stop("Sample ", input.sampleID, "not found!")
  }
  sample.qc <- subset(qc.results, qc.results$sampleID == input.sampleID)

  # Figure out the ordering of the columns in ascending order 

  if (ordering == "normal") {
    chrom.order <- c(1:22) 
  } else if(ordering == "gc") {
    # Hard coded the chromosome order by increasing GC content 
    chrom.order <- c(4,13,5,6,3,18,8,2,7,12,21,14,9,11,10,1,15,20,16,17,22,19)
  }
  pos.order <- c() 
  normals.pos.order <- c() 
  for ( i in chrom.order ) { 
    .chrname <- paste("chr", i, sep="")
    this.pos <- which(colnames(sample.qc) == .chrname)
    pos.order <- c(pos.order, this.pos)
    this.pos <- which(names(normals.ratios.mean) == .chrname)
    normals.pos.order <- c(normals.pos.order, this.pos)
  }
  sample <- as.numeric(sample.qc[1,pos.order])
  normals.mean <- as.numeric(normals.ratios.mean[normals.pos.order])
  normals.sd <- as.numeric(normals.ratios.sd[normals.pos.order])

  zscores <- (sample - normals.mean) / normals.sd
  par(mar=c(5, 5, 3, 3), pin = c(5,2.5))
  ymax <- max(abs(zscores))
  plot(1:22, zscores, pch = 3, col = "black", ylim = c(-ymax-3,ymax+3), xlab = "",xaxt = 'n') 
  affected.chrom <- which(chrom.order %in% c(13,18,21))
  points(affected.chrom, zscores[affected.chrom], pch = 17, cex = 2,col="red")
  axis(side = 1, at = 1:22, labels = names(normals.ratios.mean[normals.pos.order]), las=2)
  title(main = paste("Sample ", input.sampleID), xlab = "")

  abline(h = 3, lty = 2, col = "black", xaxt = 'n')
  abline(h = -3, lty = 2, col = "black", xaxt = 'n')

}

#############################################################
## Internal functions - not exported 
#############################################################

# Takes in a data.frame of cleaned.counts, which has length 
# of the number of chromosomes and header of the chromosome names 
# in the form of chr1, chr2, etc. 
# Returns a data.frame with the z-scores for chr21, chr18, chr13 
# the sampleIDs, callT21, callT18, callT13 and callSex 
callUnknowns <- function( cleaned.counts, sampleIDs, baselines ) { 

  message("Calling the unknowns ")
  # Load the baseline values 
  method = baselines[["method"]]
  mean21 = baselines[["mean21"]]
  mean18 = baselines[["mean18"]]
  mean13 = baselines[["mean13"]]
  meanX_females  = baselines[["meanX_females"]]
  meanY_females  = baselines[["meanY_females"]]
  meanX_males  = baselines[["meanX_males"]]
  meanY_males  = baselines[["meanY_males"]]
  
  sd21 = baselines[["sd21"]]
  sd18 = baselines[["sd18"]]
  sd13 = baselines[["sd13"]]
  sdX_females  = baselines[["sdX_females"]]
  sdY_females  = baselines[["sdY_females"]]
  sdX_males  = baselines[["sdX_males"]]
  sdY_males  = baselines[["sdY_males"]]
  
  auto.names <- rep(NULL, 22)
  for ( i in c(1:22)) { 
    auto.names[i] <- paste("chr", i, sep="")
  }
  auto.pos <- which(names(cleaned.counts) %in% auto.names)
  auto.total <- rowSums(cleaned.counts[,auto.pos])
  
  cleaned.counts$autoTotal<-auto.total
  
  if(method %in% c("zscore", "MAD")) {
    message("Using ", method)
    cleaned.counts$ratio21<- cleaned.counts$chr21/cleaned.counts$autoTotal
    cleaned.counts$ratio18<- cleaned.counts$chr18/cleaned.counts$autoTotal
    cleaned.counts$ratio13<- cleaned.counts$chr13/cleaned.counts$autoTotal
    cleaned.counts$ratioX<- cleaned.counts$chrX/cleaned.counts$autoTotal
    cleaned.counts$ratioY<- cleaned.counts$chrY/cleaned.counts$autoTotal
  } else if(method == "NCV") {
    message("Using NCV")
    cleaned.counts$ratio21<- cleaned.counts$chr21/cleaned.counts$autoTotal
    cleaned.counts$ratio18<- cleaned.counts$chr18/cleaned.counts$chr8
    cleaned.counts$ratio13<- cleaned.counts$chr13/(cleaned.counts$chr4 + cleaned.counts$chr5)
    cleaned.counts$ratioX <- cleaned.counts$chrX/(cleaned.counts$chr3 + cleaned.counts$chr4)
    cleaned.counts$ratioY<- cleaned.counts$chrY/cleaned.counts$autoTotal    
  } 
  
  cleaned.counts$z21<- ( cleaned.counts$ratio21 - mean21 ) / sd21
  cleaned.counts$z18<- ( cleaned.counts$ratio18 - mean18 ) / sd18
  cleaned.counts$z13<- ( cleaned.counts$ratio13 - mean13 ) / sd13 
  cleaned.counts$zX_females<- ( cleaned.counts$ratioX - meanX_females ) / sdX_females
  cleaned.counts$zY_females<- ( cleaned.counts$ratioY - meanY_females ) / sdY_females
  cleaned.counts$zX_males<- ( cleaned.counts$ratioX - meanX_males ) / sdX_males
  cleaned.counts$zY_males<- ( cleaned.counts$ratioY - meanY_males ) / sdY_males 
  
  cleaned.counts$callT21 <- ( cleaned.counts$z21 > 3 )
  cleaned.counts$callT18 <- ( cleaned.counts$z18 > 3 )
  cleaned.counts$callT13 <- ( cleaned.counts$z13 > 3 )
  
  cleaned.counts$callSex <- cleaned.counts$callT21 
  
  CV_chrX <- sdX_females / meanX_females 
 
  for ( i in 1:nrow(cleaned.counts)) {
    if (cleaned.counts[i,"zY_females"] > 3) {
       cleaned.counts[i, "callSex"] <- "Male"
       cleaned.counts[i, "FF_chrX"] <- abs(2*cleaned.counts[i, "zX_females"] * CV_chrX)  
       cleaned.counts[i, "ratio_chrY"] <- cleaned.counts[i,"chrY"]/cleaned.counts[i, "autoTotal"]
    } else if (cleaned.counts[i,"zX_females"] > -2 & cleaned.counts[i, "zY_females"] < 3 ) {
      cleaned.counts[i,"callSex"] <- "Female"
    } else if (cleaned.counts[i,"zX_females"] < -3 & cleaned.counts[i, "zY_females"] < 3) {
      cleaned.counts[i,"callSex"] <- "Turner"
    } else {
       cleaned.counts[i, "callSex"] <- "No call"
    }
  }
   
  # For samples called as male, also provide a likelihood ratio using the 
  # estimated fetal fraction 
  males <- which(cleaned.counts$callSex == "Male") 
  for (m in males) {
    ratio21 <- cleaned.counts$ratio21[m]
    ff      <- cleaned.counts$FF_chrX[m]
    likeli_notri <- dnorm(ratio21, mean = mean21, sd = sd21)
    
    ff_rand <- rnorm(10000, mean = ff, sd = 0.02)
    likeli_tri_list <- seq(1, 10000)
    for (i in 1:10000) {
      likeli_tri_list[i] <- dnorm(ratio21, mean = mean21*(1+0.5*ff_rand[i]), sd = sd21)
    }
    likeli_tri   <- mean(likeli_tri_list) 
    
    likeli_ratio <- log10(likeli_notri/likeli_tri)
    cleaned.counts[m, "log_li"] <- likeli_ratio 
  }

  # If there are no males called, still add the log_li column 
  if (length(males) == 0) { 
    cleaned.counts[,"log_li"] <- rep(NA, nrow(cleaned.counts)) 
    cleaned.counts[,"FF_chrX"] <- rep(NA, nrow(cleaned.counts)) 
    cleaned.counts[,"ratio_chrY"] <- rep(NA, nrow(cleaned.counts)) 
  } 

  results.cols <- c("z21", "z18", "z13", "zX_females", "zY_males", "zY_females", 
                    "FF_chrX", "ratio_chrY", "callT21", "callT18", "callT13",
                    "callSex", "log_li")
  results.df <- cleaned.counts[,results.cols]  
  results.df$SampleIDs <- sampleIDs 
  results.df <- results.df[,c(14, 1:13)]
 
  return(results.df)  
}

# Returns a data.frame with the ratio of counts mapped for every chromosome, 
# the total read counts, and whether QC was passed 
data.qc <- function(cleaned.counts, sampleIDs, baselines) {
  normals.ratios.mean <- baselines[["normals.mean"]]
  normals.ratios.sd   <- baselines[["normals.sd"]]    

  auto.names <- rep(NULL, 22)
  j <- 1 
  for ( i in c(1:22) ) { 
    auto.names[j] <- paste("chr", i, sep="")
    j <- j + 1 
  }

  auto.pos <- which(names(cleaned.counts) %in% auto.names)
  auto.total <- rowSums(cleaned.counts[,auto.pos])
  cleaned.counts$autoTotal<-auto.total  
  cleaned.ratios <- cleaned.counts[,auto.pos] / cleaned.counts$autoTotal 

  qc.total.counts <- data.frame(QC.reads = (cleaned.counts$autoTotal > 2e6))
  qc.total.counts["totalAutoCounts"] <- cleaned.counts$autoTotal 

  qc.ref.match <- c() 
  for (i in 1:nrow(cleaned.ratios)) {

    cleaned.ratios.zscores <- (cleaned.ratios[i,] - normals.ratios.mean) / normals.ratios.sd
    .qc.ref.match <- data.frame(cleaned.ratios[i,]) 
    .qc.ref.match["sampleID"] <- sampleIDs[i]
    chr.excl <- which(names(cleaned.ratios.zscores) %in% c("chr18", "chr13", "chr21"))
    flag.chr <- which(abs(cleaned.ratios.zscores[-chr.excl]) > 3) 

    if (length(flag.chr) > 0) {
      message("Sample ", sampleIDs[i], " has ", length(flag.chr), 
              " chromosome with counts outside the expected range")
      .qc.ref.match["QC.refmatch"] <- FALSE 
    } else {
      .qc.ref.match["QC.refmatch"] <- TRUE 
    }

    qc.ref.match <- rbind(qc.ref.match, .qc.ref.match) 
  }

  qc.results <- cbind(qc.ref.match, qc.total.counts)
  
  return(qc.results)
}
