#' @title createReferenceSetFromCounts
#' 
#' @description
#' This function creates a reference set from a binned counts file  
#'  
#' @param binned.counts.file file name of the binned counts. The binned counts file 
#'        should be comma delimited, and the first line need to be the chromosome 
#'        names of each bin 
#' @param outcomes data.frame with column names: Dx, Gender, SampleID 
#' @param combined.counts.fname file name to write to for the combined counts
#'        per chromosomes, default is not to write result to file  
#' @param method either "zscore", "NCV" or "MAD", default is zscore         
#' @param gcCorrect whether to do gc correction or not (True = do the correction)
#' @param PCA whether to do PCA correction or not (True = do the correction)
#' @param numPC number of principal components to discard (default = 10)  
#' @param masked.counts.file file name of the masked counts file 
#' @param gcContentFile file name of a Rdata object with the gcContent data 
#' @param filterBin whether to filter bins based on unusually high counts and high
#'        variance, default is to filter 
#' @param removeOutlierSamples whether to remove samples which has a low correlation value 
#'        to the rest of the reference set, default is FALSE 
#' @param cleaned.binned.counts.fname file name to write to for the corrected binned
#'        counts, default is not to write to file 
#'
#' @return class of rapidr.ref which can then be used to test unknown samples 
#' @export 
#' @importFrom data.table fread  
#'  
createReferenceSetFromCounts <- function( binned.counts.file, outcomes,
                                          combined.counts.fname = NULL, 
                                        method = "zscore", gcCorrect = FALSE, 
                                          gcContentFile = NULL, filterBin = TRUE,
                                          removeOutlierSamples = FALSE, 
                                          PCA = FALSE, numPC = 10, 
                                          masked.counts.file = NULL,
                                          cleaned.binned.counts.fname = NULL) {


    ptm <- proc.time() 
    #print(ptm) 
    message("Loading binned counts file")
    masked.binned.counts <- NULL 

    # Make sure there is a gcContentFile provided if doing gc correction 
    if (gcCorrect == TRUE) {
      if (is.null(gcContentFile)) {
        stop("You need to provide a gcContent file. Try running makeGCContentData()")
      }
    }
    
    binned.counts <- fread(binned.counts.file, sep=',', colClasses=list(character=1), header=TRUE)
    
    # TO FIX: This is not the best way to subset the data.table
    # But the way that is advised in the help pages seem to give a strange bug 
    #n.cols <- ncol(binned.counts)
    #ind<-rep(FALSE, n.cols)
    #ind[1] <-TRUE

    bin.names <- names(binned.counts)
    bin.names <- bin.names[2:length(bin.names)] 

    #message("making data.frame") 
    binned.counts <- as.data.frame(binned.counts) 
    sampleIDs <- binned.counts[,1] 
    binned.counts <- binned.counts[,-c(1)]  

    if (!is.null(masked.counts.file)) {
      masked.binned.counts <- fread(masked.counts.file, sep=',', colClasses=list(character=1))
      masked.binned.counts <- data.frame(masked.binned.counts) 
      masked.binned.counts <- masked.binned.counts[,-c(1)] 
    }
   
    n.all.samples <- nrow(binned.counts)
    n.bins <- ncol(binned.counts)
    binned.counts <- data.matrix(binned.counts) 

    if (!is.null(masked.counts.file)) {
       masked.binned.counts <- as.matrix(masked.binned.counts)  
    }
    
    #message("Memory used: ") 
    #print(sort(sapply(ls(),function(x){object.size(get(x))}) ))  

    if ( gcCorrect & !is.null(masked.counts.file)) {
      message("Doing gc correction")
      binned.counts <- gcCorrectCounts(binned.counts, bin.names, masked.binned.counts, gcContentFile)
    } else if ( gcCorrect & is.null(masked.counts.file)) {
      binned.counts <- gcCorrectCounts(binned.counts, bin.names, binned.counts, gcContentFile)
    } else if (!is.null(masked.counts.file)) {
      # If there is a masked binned counts file, use that 
      binned.counts <- masked.binned.counts 
    }
    
    discard.pos <- c()
    n.discard <- 0  
    
    #message("Memory used: ") 
    #print(sort(sapply(ls(),function(x){object.size(get(x))}) ))  

    # Checking every sampleID has an outcome 
    message("Checking every sampleID has an outcome")
    for (i in 1:n.all.samples) { 
      eachSample <- sampleIDs[i]
      pos <- which(outcomes$SampleID == as.character(eachSample))
      if (length(pos) == 0) { 
        message("No outcomes for Sample ", eachSample)    
        n.discard <- n.discard + 1 
        discard.pos[n.discard] <- i 
      }
    }
    
    # Match outcomes data to the sampleIDs 
    sampleIDs.with.outcomes <- data.frame(sampleIDs = sampleIDs)
    for ( i in 1:length(sampleIDs)) { 
      pos <- which(outcomes$SampleID == as.character(sampleIDs[i]))    
      if (length(pos) == 0) { 
        next 
      }
      sampleIDs.with.outcomes[i,"Dx"] <- outcomes[pos, "Dx"]
      sampleIDs.with.outcomes[i,"Gender"] <- outcomes[pos, "Gender"]
    }
   
    # Check we have at least two female and two male samples 
    n.females <- length(which(sampleIDs.with.outcomes[,"Gender"] == "Female") ) 
    n.males <- length(which(sampleIDs.with.outcomes[,"Gender"] == "Male") ) 
    
    if (n.females < 3 | n.males < 3) { 
      stop("Requires at least 2 samples of females and 2 samples of males.")        
    }  

    # Find bins to exclude in chromosome Y 
    message("Finding bins to exclude in Chr Y ")
    if(length(discard.pos) > 0) {
      excl.bins <- find.chrY.excl.bins(binned.counts, bin.names, sampleIDs.with.outcomes[-discard.pos,]) 
    } else {
      excl.bins <- find.chrY.excl.bins(binned.counts, bin.names, sampleIDs.with.outcomes) 
    }

    total.counts <-rowSums(binned.counts) 
    binned.ratios <- binned.counts 
    
    for (i in 1:n.all.samples) { 
      this.total<-total.counts[i] 
      if (this.total < 2e6) { 
         message("Sample ", sampleIDs[i], " has less than 2 million counts")
         n.discard <- n.discard + 1 
         discard.pos[n.discard] <- i 
         next 
      }
      # Multiply by 1million to avoid small ratio numbers 
      binned.ratios[i,] <- binned.ratios[i,]/this.total * 1e6     
    }
       
    message("Number of discarded samples: ", n.discard)
    normals.ids <- outcomes[which(outcomes$Dx == "Normal"), "SampleID"] 
    normals.pos <- which(sampleIDs %in% normals.ids)
    normals.pos <- subset(normals.pos, ! normals.pos %in% discard.pos) 

    females.ids <- outcomes[which(outcomes$Dx == "Normal" & outcomes$Gender == "Female"), "SampleID"] 
    females.pos <- which(sampleIDs %in% females.ids)
    females.pos <- subset(females.pos, ! females.pos %in% discard.pos) 
    
    # Find bins to exclude in other chromosomes 
    if(filterBin) { 
      message("Finding bins to exclude in other chromosomes")
      if(length(discard.pos) > 0) {
        other.excl.bins <- find.excl.bins(binned.ratios[normals.pos,], bin.names)
      } else {
        other.excl.bins <- find.excl.bins(binned.ratios[normals.pos,], bin.names) 
      }
      excl.bins <- c(excl.bins, other.excl.bins) 
    }   

    binned.counts[,excl.bins] <- 0.0 

    # Re-calculate the bin ratios to take into account the bins that has been excluded 
    total.counts <-rowSums(binned.counts) 
    binned.ratios <- binned.counts 
    for (i in 1:n.all.samples) { 
      this.total<-total.counts[i] 
      binned.ratios[i,] <- binned.ratios[i,]/this.total * 1e6     
    }

    rm(binned.counts)
    rm(masked.binned.counts)

    if (length(discard.pos) > 0) {
      good.ratios    <- binned.ratios[ -discard.pos, ]
      good.sampleIDs <- sampleIDs[-discard.pos]
    } else {
      good.ratios    <- binned.ratios
      good.sampleIDs <- sampleIDs      
    }
       
    if ( PCA ) {
      #data.mat.clean.for.PCA <- binned.ratios[ normals.pos, ]
      #PCA_output <- doPCA(binned.ratios[normals.pos, ]) 
      if (numPC < 1) { 
         stop("Number of PC to remove must be at least 1 up to the number of female samples") 
      } 
      if (length(females.pos) < numPC ) { 
         stop("There are less female samples than the number of principal components. PCA error correction will not work properly. Try reducing the numPC or increase the number of samples") 
      } 

      data.mat.clean.for.PCA <- binned.ratios[ females.pos, ]
      PCA_output <- doPCA(binned.ratios[females.pos, ], numPC = numPC) 
      #save(PCA_output, file = "~/UCL/PhaseI_all/PCA_output.Rdata")
      #load("~/UCL/PhaseI_all/PCA_output.Rdata")
      # Use the PCA results to correct the counts 
      if(length(discard.pos) > 0) {
        cleaned.good.ratios <- correct.counts.with.PCA(PCA_output, binned.ratios[-discard.pos, ])
      } else {
        cleaned.good.ratios <- correct.counts.with.PCA(PCA_output, binned.ratios)
      }
    } else {
      # Don't do any corrections 
      message("Not doing PCA")
      if (length(discard.pos) > 0) {
         cleaned.good.ratios <- binned.ratios[-discard.pos,]
      } else {
        cleaned.good.ratios <- binned.ratios        
      }
    }
      
    # Find the mean of each bin from the normals 
    if (length(discard.pos) > 0) {
      cleaned.normals.pos <- which(sampleIDs[-discard.pos] %in% normals.ids)
    } else {
      cleaned.normals.pos <- normals.pos 
    }
    bin.means <- colMeans(cleaned.good.ratios[cleaned.normals.pos,])
    bin.sd    <- apply(cleaned.good.ratios[cleaned.normals.pos,], 2, sd)

    # Checking for samples which are outliers 
    message("Checking for outlier samples...") 
    ratios.cor <- cor(bin.means, t(cleaned.good.ratios)) 
    if (length(discard.pos) > 0) { 
       ratios.df <- data.frame(corr = ratios.cor[1,], sampleID = sampleIDs[-discard.pos]) 
    } else { 
       ratios.df <- data.frame(corr = ratios.cor[1,], sampleID = sampleIDs) 
    } 
    ratios.order <- order(ratios.df[,1]) 
    low.corr.samples <- ratios.df[which(ratios.df$corr < 0.8), 2] 
    
    if (removeOutlierSamples == FALSE) { 
       message("Samples with correlation to the mean bin counts of < 0.8. You may want to remove these from the reference set.") 
       print(low.corr.samples) 
    } else {
       if(length(discard.pos) >0) {
          cleaned.good.ratios <- cleaned.good.ratios[-which(sampleIDs[-discard.pos] %in% low.corr.samples),]
       } else { 
          cleaned.good.ratios <- cleaned.good.ratios[-which(sampleIDs %in% low.corr.samples),]
       } 
       new.discards <- which(sampleIDs %in% low.corr.samples) 
       discard.pos <- c(discard.pos, new.discards) 
       n.outliers <- length(new.discards) 
       message("... removing ", n.outliers, " samples as outliers") 

       # Re-calculate the bin.means and bin.sd after taking out the discards  
       cleaned.normals.pos <- which(sampleIDs[-discard.pos] %in% normals.ids) 
       bin.means <- colMeans(cleaned.good.ratios[cleaned.normals.pos,])
       bin.sd    <- apply(cleaned.good.ratios[cleaned.normals.pos,], 2, sd)
       
    } 

    message("Summing counts per chromosome...")
    if (length(discard.pos) > 0) {
      cleaned.counts.per.chr <- sum.counts(cleaned.good.ratios, bin.names, total.counts[-discard.pos], 
                                         sampleIDs[-discard.pos])   
      baselines <- findBaseline(cleaned.counts.per.chr, sampleIDs.with.outcomes[-discard.pos,], method = method)
    } else {
      cleaned.counts.per.chr <- sum.counts(cleaned.good.ratios, bin.names, total.counts, 
                                           sampleIDs)      
      baselines <- findBaseline(cleaned.counts.per.chr, sampleIDs.with.outcomes, method = method)      
    }
    
    
    # If there is a file name provided, write the cleaned, combined counts to a file 
    if ( ! is.null(combined.counts.fname) ) {
      write.csv(cleaned.counts.per.chr, file = combined.counts.fname, quote = FALSE, row.names = FALSE)     
    }
    
    # Put the relevant data of the reference set somewhere 
    ref.data.set <- list()  
    
    if ( PCA ) {
      ref.data.set[["PCA"]] <- PCA_output      
    } else { 
      ref.data.set[["PCA"]] <- NULL 
    }

    ref.data.set[["baselines"]] <- baselines 
    ref.data.set[["excl.bins"]] <- excl.bins 
    ref.data.set[["do.gcCorrect"]] <- gcCorrect 
    ref.data.set[["do.PCA"]] <- PCA 
    ref.data.set[["bin.means"]] <- bin.means
    ref.data.set[["bin.sd"]] <- bin.sd
    
    class(ref.data.set) <- "rapidr.ref"

    # Evaluate the reference set 
    refset.calls <- callUnknowns(cleaned.counts.per.chr, cleaned.counts.per.chr$SampleID, baselines)
    refset.results <- evalPerformance(refset.calls, outcomes)
    ref.data.set[["baseline.perf"]] <- refset.results 
    
    # Write the cleaned binned counts to a file if there is a file name 
    # provided 
    if ( ! is.null(cleaned.binned.counts.fname) ) {
      message("Writing the cleaned binned counts to the file provided.")
      # Convert ratios back to counts 
      if (length(discard.pos) > 0) { 
         cleaned.total <- total.counts[-discard.pos]
      } else { 
         cleaned.total <- total.counts
      } 
      binned.counts <- cleaned.good.ratios
      for (i in 1:nrow(cleaned.good.ratios)) { 
        this.total<-cleaned.total[i]
        # Multiply by 1million to avoid small ratio numbers 
        binned.counts[i,] <- cleaned.good.ratios[i,] * this.total / 1e6     
      }
      rm(cleaned.good.ratios) 
      sampleIDs <- as.character(sampleIDs)
      if (length(discard.pos) > 0 ) {  
         binned.counts <- cbind(sampleIDs[-discard.pos], binned.counts)
      } else { 
         binned.counts <- cbind(sampleIDs, binned.counts)
      } 
      write.table(binned.counts, file = cleaned.binned.counts.fname, sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE)     
    }    
   
    #message("Time used:")  
    #print(proc.time() - ptm) 
    return(ref.data.set)
    
}

#' @title writeCleanedCountsFiles
#' 
#' @description
#' This function takes a binned.counts.file and applies GC correction or PCA correction 
#' and writes out the resulting binned counts as a a new file   
#'  
#' @param binned.counts.file file name of the binned counts. The binned counts file 
#'        should be comma delimited, and the first line need to be the chromosome 
#'        names of each bin       
#' @param gcContentFile file name of a Rdata object with the gcContent data 
#' @param gcCorrect whether to do gc correction or not (True = do the correction)
#' @param PCA whether to do PCA correction or not (True = do the correction)
#' @param cleaned.binned.counts.fname file name to write to for the corrected binned
#'        counts, default is not to write to file 
#' 
#' @export 
#' @importFrom data.table fread 
#'  
writeCleanedCountsFile <- function( binned.counts.file, cleaned.binned.counts.fname, 
                                    gcContentFile, 
                                          gcCorrect = FALSE, 
                                          PCA = FALSE) {
  
  message("Loading binned counts file")
  
  # Make sure there is a gcContentFile provided if doing gc correction 
  if (gcCorrect == TRUE) {
    if (is.null(gcContentFile)) {
      stop("You need to provide a gcContent file. Try running makeGCContentData()")
    }
  }
  
  binned.counts <- fread(binned.counts.file, sep=',', colClasses=list(character=1), header=TRUE)
  
  # TO FIX: This is not the best way to subset the data.table
  # But the way that is advised in the help pages seem to give a strange bug 
  n.cols <- ncol(binned.counts)
  #ind<-rep(FALSE, n.cols)
  #ind[1] <-TRUE
  
  binned.counts <- data.frame(binned.counts) 
  sampleIDs <- binned.counts[,1] 
  binned.counts <- binned.counts[,-c(1)]     

  #sampleIDs<-data.frame(binned.counts[,ind,with=FALSE])
  #sampleIDs <- sampleIDs[,1] 
  #binned.counts <- binned.counts[,!ind, with=FALSE]
  bin.names <- names(binned.counts)

  n.all.samples <- nrow(binned.counts)
  n.bins <- ncol(binned.counts)
  
  binned.counts <- as.matrix(binned.counts) 
  masked.binned.counts <- as.matrix(binned.counts) 
  
  if ( gcCorrect ) {
    message("Doing gc correction")
    binned.counts <- gcCorrectCounts(binned.counts, bin.names, masked.binned.counts, gcContentFile)
  } else  {
    binned.counts <- gcCorrectCounts(binned.counts, bin.names, binned.counts, gcContentFile)
  } 
  
  discard.pos <- c()
  n.discard <- 0  
  
  total.counts <-rowSums(binned.counts) 
  binned.ratios <- binned.counts 
  rm(binned.counts)
  
  for (i in 1:n.all.samples) { 
    this.total<-total.counts[i] 
    if (this.total < 2e6) { 
      message("Sample ", sampleIDs[i], " has less than 2 million counts")
      n.discard <- n.discard + 1 
      discard.pos[n.discard] <- i 
      next 
    }
    # Multiply by 1million to avoid small ratio numbers 
    binned.ratios[i,] <- binned.ratios[i,]/this.total * 1e6     
  }
  
  message("Number of discarded samples: ", n.discard) 
  
  if (length(discard.pos) > 0) {
    good.ratios    <- binned.ratios[ -discard.pos, ]
    good.sampleIDs <- sampleIDs[-discard.pos]
  } else {
    good.ratios    <- binned.ratios
    good.sampleIDs <- sampleIDs      
  }
  
  if ( PCA ) {
    data.mat.clean.for.PCA <- binned.ratios
    PCA_output <- doPCA(binned.ratios) 
    # Use the PCA results to correct the counts 
    if(length(discard.pos) > 0) {
      cleaned.good.ratios <- correct.counts.with.PCA(PCA_output, binned.ratios[-discard.pos, ])
    } else {
      cleaned.good.ratios <- correct.counts.with.PCA(PCA_output, binned.ratios)
    }
    gc()
  } else {
    # Don't do any corrections 
    message("Not doing PCA")
    if (length(discard.pos) > 0) {
      cleaned.good.ratios <- binned.ratios[-discard.pos,]
    } else {
      cleaned.good.ratios <- binned.ratios        
    }
  }
  
  # Write the cleaned binned counts to a file if there is a file name 
  # provided 
  if ( ! is.null(cleaned.binned.counts.fname) ) {
    message("Writing the cleaned binned counts to the file provided.")
    # Convert ratios back to counts 
    if ( length(discard.pos) > 0 ) {
      cleaned.total <- total.counts[-discard.pos]
    } else {
      cleaned.total <- total.counts
    }
    binned.counts <- cleaned.good.ratios
    for (i in 1:nrow(cleaned.good.ratios)) { 
      this.total<-cleaned.total[i]
      # Multiply by 1million to avoid small ratio numbers 
      binned.counts[i,] <- cleaned.good.ratios[i,] * this.total / 1e6     
    }
    rm(cleaned.good.ratios) 
    sampleIDs <- as.character(sampleIDs)
    if ( length(discard.pos ) > 0 ) {
       binned.counts <- cbind(sampleIDs[-discard.pos], binned.counts)
    } else {
      binned.counts <- cbind(sampleIDs, binned.counts)      
    }
    write.table(binned.counts, file = cleaned.binned.counts.fname, sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE)     
  }    
  
}

#' @title makeBinnedCountsFile
#' 
#' @description
#' This function takes in a list of bam files and creates a binned   
#' counts file. If a mask file is provided, it will also create 
#' a masked binned counts file. The output file is comma separated and 
#' the first column is the sampleID, the header is the chromosome name of 
#' each bin. 
#' 
#' @param bam.file.list list of bam file names 
#' @param sampleIDs list of sampleIDs, assumed to be in the same order as 
#'        the bam files in bam.file.list 
#' @param binned.counts.fname file name of the output binned counts file 
#' @param mask  file name of a bed file with the regions to be masked out. Default
#'        is no mask file 
#' @param k bin size in number of bases. Default is 20,000 bp 
#' 
#' 
#' @export 
#' @seealso \code{\link{createReferenceSetFromCounts}}
#' 

makeBinnedCountsFile <- function (bam.file.list, sampleIDs, binned.counts.fname, mask = NULL, k = 20000) {
  message("Binning counts in bam files")
  res <- BinListOfBam(bam.file.list, mask = mask, k = k)
  all.binnedCounts <- res[[1]]
  masked.binnedCounts <- res[[2]]
  #write.table(colnames(all.binnedCounts), file = bin.names.fname, quote = FALSE, row.names = FALSE, col.names = FALSE)
  bin.names <- colnames(all.binnedCounts)
  bin.names <- t(bin.names)
  rownames(bin.names) <- c("SampleID")
  write.table(bin.names, file = binned.counts.fname, quote = FALSE, row.names = TRUE, col.names=FALSE, sep=',')  
  all.binnedCounts <- as.matrix(all.binnedCounts)
  rownames(all.binnedCounts) <- sampleIDs
  write.table(all.binnedCounts, file = binned.counts.fname, append = TRUE, quote = FALSE, row.names = TRUE, col.names=FALSE, sep=',')  

  masked.counts.fname <- paste(binned.counts.fname, ".mask", sep="") 
  masked.binnedCounts <- as.matrix(masked.binnedCounts)
  write.table(bin.names, file = masked.counts.fname, quote = FALSE, row.names = TRUE, col.names=FALSE, sep=',')  
  rownames(masked.binnedCounts) <- sampleIDs
  #write.table(masked.binnedCounts, file = masked.counts.fname, append = TRUE, quote = FALSE, row.names = TRUE, col.names=FALSE, sep=',')  
}


#############################################################
## Internal functions - not exported 
#############################################################

# Internal function for summing the counts for each chromosome given the 
# counts ratios 
sum.counts <- function ( cleaned.good.ratios, bin.names, total.counts, sampleIDs ) { 

    # Check that number of samples in the ratios matrix is the same as in total.counts 
    if (nrow(cleaned.good.ratios) != length(total.counts) ) {
      stop("Number of samples in ratios matrix not the same as in total.counts? Please check.")
    }

    # Check that number of bins in the ratios matrix is the same as in bin.names
    if (ncol( cleaned.good.ratios) != length(bin.names) ) {
      stop("Number of bins ratios matrix not the same as in bin.names? Please check.")
    }

    chr.names <- levels(as.factor(bin.names) )
    #print(length(chr.names))
    
    counts.per.chr <- data.frame(matrix(ncol = (length(chr.names) + 1), nrow = nrow(cleaned.good.ratios)))

    counts.per.chr[,1] <- sampleIDs
    names(counts.per.chr)[1] <- "SampleID"
    i <- 2 
    for ( each.chr in chr.names ) {
      chr.pos <- which(bin.names == each.chr)
      for ( j in 1:nrow(cleaned.good.ratios)) {
        counts.per.chr[j,i] <- sum(cleaned.good.ratios[j,chr.pos]) * total.counts[j] / 1e6       
      }

      names(counts.per.chr)[i] <- each.chr
      i <- i + 1 
    }

    return(counts.per.chr)
}

findBaseline <- function(cleaned.counts, sampleIDs.with.outcomes, method = "zscore" ) {

   message("Finding baseline values...")
   if (nrow(cleaned.counts) != nrow(sampleIDs.with.outcomes) ) {
     stop("Ratios matrix has a different number of samples to the outcomes table? Check.")
   }
   
   if (!method %in% c("zscore", "NCV", "MAD")) {
     stop("Method needs to be one of zscore, NCV or MAD")
   }
   
   normals.ids <- sampleIDs.with.outcomes[which(sampleIDs.with.outcomes$Dx == "Normal"), "SampleID"] 
   normals <- cleaned.counts[which(sampleIDs.with.outcomes$Dx == "Normal"), ]

   females.ids <- sampleIDs.with.outcomes[which(sampleIDs.with.outcomes$Dx == "Normal" & 
                                                  sampleIDs.with.outcomes$Gender == "Female"), "SampleID"] 
   females <- cleaned.counts[which(sampleIDs.with.outcomes$Dx == "Normal" & 
                                     sampleIDs.with.outcomes$Gender == "Female"), ]

   males.ids <- sampleIDs.with.outcomes[which(sampleIDs.with.outcomes$Dx == "Normal" & 
                                                  sampleIDs.with.outcomes$Gender == "Male"), "SampleID"] 
   males <- cleaned.counts[which(sampleIDs.with.outcomes$Dx == "Normal" & 
                                     sampleIDs.with.outcomes$Gender == "Male"), ]
   
   auto.names <- rep(NULL, 22)
   for ( i in c(1:22) ) { 
      auto.names[i] <- paste("chr", i, sep="")
   }
   auto.pos <- which(names(cleaned.counts) %in% auto.names)
   auto.total <- rowSums(normals[,auto.pos])

   # Subset of the autosomes which does not include 
   # chromosome with possible trisomies (chr21, chr18, ch13)
   sub.auto.names <- rep(NULL, 19)
   j <- 1 
   for ( i in c(1:12, 14:17, 19:20, 22) ) { 
     sub.auto.names[j] <- paste("chr", i, sep="")
     j <- j + 1 
   }
   sub.auto.pos <- which(names(cleaned.counts) %in% sub.auto.names)
   sub.auto.total <- rowSums(normals[,sub.auto.pos])

   auto.total.females <- rowSums(females[,auto.pos])
   auto.total.males <- rowSums(males[,auto.pos])
   
   normals$autoTotal<-auto.total
   normals$sub.autoTotal <- sub.auto.total
   females$autoTotal <- auto.total.females 
   males$autoTotal <- auto.total.males    
   
   normals.ratios <- normals[,auto.pos] / normals$autoTotal
   normals.ratios.mean <- apply(normals.ratios, 2, mean)
   normals.ratios.sd <- apply(normals.ratios, 2, sd)
  
   
   if (method == "zscore" ) {
     normals$ratio21<- normals$chr21/normals$autoTotal
     normals$ratio18<- normals$chr18/normals$autoTotal
     normals$ratio13<- normals$chr13/normals$autoTotal
     normals$ratioX <- normals$chrX/normals$autoTotal 
     normals$ratioY <- normals$chrY/normals$autoTotal
     
     females$ratioX <- females$chrX/females$autoTotal 
     females$ratioY <- females$chrY/females$autoTotal 
     males$ratioX <- males$chrX/males$autoTotal 
     males$ratioY <- males$chrY/males$autoTotal
     
     mean21<-mean(normals$ratio21)
     mean18<-mean(normals$ratio18)
     mean13<-mean(normals$ratio13)
   
     meanX_females <- mean(females$ratioX)
     meanY_females <- mean(females$ratioY)

     meanX_males <- mean(males$ratioX)
     meanY_males <- mean(males$ratioY)
     
     sd21 <- sd(normals$ratio21)
     sd18 <- sd(normals$ratio18)
     sd13 <- sd(normals$ratio13)
   
     sdX_females <- sd(females$ratioX) 
     sdY_females <- sd(females$ratioY) 
     sdX_males <- sd(males$ratioX) 
     sdY_males <- sd(males$ratioY) 
     
     baselines <- list() 
     baselines[["method"]] <- "zscore"
     
   } else if (method == "NCV") { 
     normals$ratio21<- normals$chr21/normals$autoTotal
     normals$ratio18<- normals$chr18/normals$chr8
     normals$ratio13<- normals$chr13/(normals$chr4 + normals$chr5)
     normals$ratioX <- normals$chrX/(normals$chr3 + normals$chr4) 
     normals$ratioY <- normals$chrY/normals$autoTotal 
     
     females$ratioX <- females$chrX/(females$chr3 + females$chr4) 
     females$ratioY <- females$chrY/females$autoTotal 

     males$ratioX <- males$chrX/(males$chr3 + males$chr4) 
     males$ratioY <- males$chrY/males$autoTotal 
     
     mean21<-mean(normals$ratio21)
     mean18<-mean(normals$ratio18)
     mean13<-mean(normals$ratio13)
     
     meanX_females <- mean(females$ratioX)
     meanY_females <- mean(females$ratioY)

     meanX_males <- mean(males$ratioX)
     meanY_males <- mean(males$ratioY)
     
     sd21 <- sd(normals$ratio21)
     sd18 <- sd(normals$ratio18)
     sd13 <- sd(normals$ratio13)
     
     sdX_females <- sd(females$ratioX) 
     sdY_females <- sd(females$ratioY)  
     sdX_males <- sd(males$ratioX) 
     sdY_males <- sd(males$ratioY) 
     
     baselines <- list() 
     baselines[["method"]] <- "NCV"
   } else if (method == "MAD") {
     normals$ratio21<- normals$chr21/normals$autoTotal
     normals$ratio18<- normals$chr18/normals$autoTotal    
     normals$ratio13<- normals$chr13/normals$autoTotal
     normals$ratioX <- normals$chrX/normals$autoTotal 
     normals$ratioY <- normals$chrY/normals$autoTotal
     
     females$ratioX <- females$chrX/females$autoTotal 
     females$ratioY <- females$chrY/females$autoTotal 
     males$ratioX <- males$chrX/males$autoTotal 
     males$ratioY <- males$chrY/males$autoTotal
     
     mean21<-median(normals$ratio21)
     mean18<-median(normals$ratio18)
     mean13<-median(normals$ratio13)
     
     meanX_females <- median(females$ratioX)
     meanY_females <- median(females$ratioY)
     
     meanX_males <- median(males$ratioX)
     meanY_males <- median(males$ratioY)
     
     sd21 <- mad(normals$ratio21)
     sd18 <- mad(normals$ratio18)
     sd13 <- mad(normals$ratio13)
     
     sdX_females <- mad(females$ratioX) 
     sdY_females <- mad(females$ratioY) 
     sdX_males <- mad(males$ratioX) 
     sdY_males <- mad(males$ratioY) 
     
     baselines <- list() 
     baselines[["method"]] <- "MAD"
     
   } else {
     message("Method unknown! Need to be either zscore, NCV or MAD.")     
   }
   
   baselines[["mean21"]] <- mean21 
   baselines[["mean18"]] <- mean18 
   baselines[["mean13"]] <- mean13 
   baselines[["meanX_females"]]  <- meanX_females 
   baselines[["meanY_females"]]  <- meanY_females  
   baselines[["meanX_males"]]  <- meanX_males 
   baselines[["meanY_males"]]  <- meanY_males  
   
   baselines[["sd21"]] <- sd21 
   baselines[["sd18"]] <- sd18 
   baselines[["sd13"]] <- sd13 
   baselines[["sdX_females"]]  <- sdX_females
   baselines[["sdY_females"]]  <- sdY_females     
   baselines[["sdX_males"]]  <- sdX_males
   baselines[["sdY_males"]]  <- sdY_males 
   
   baselines[["normals.mean"]] <- normals.ratios.mean
   baselines[["normals.sd"]]   <- normals.ratios.sd 
   
   return(baselines)
}

