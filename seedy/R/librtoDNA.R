librtoDNA <-
function(sampleID, libr, nuc, ref.strain, key, sampletime=NULL, strings=FALSE, 
         filename=NULL, format="nexus") {
  ngenomes <- length(sampleID)
  DNAoutput <- NULL
  if (!is.null(filename)) {
    strings <- TRUE
  }
  ntides <- c("C", "A", "G", "T")
  if (!is.null(filename) && format=="nexus") {
    write(paste("#nexus\nbegin data;\ndimensions ntax=", length(sampleID), 
        " nchar=", length(ref.strain), ";\nformat datatype=dna symbols=\"CAGT\" missing=? gap=-;\nmatrix\n", sep=""),
        file=filename)
    strings <- TRUE
  }
  if (format=="fasta") {
    strings <- FALSE
  }
  for (i in 1:ngenomes) {
    K <- ref.strain
    if (!is.na(nuc[[which(key==sampleID[i])]][1])) {
      K[ libr[[which(key==sampleID[i])]] ] <- nuc[[which(key==sampleID[i])]]
    }
    for (j in 1:length(ntides)) {
      K[which(K==j)] <- ntides[j]
    }
    if (strings) {
      DNAoutput <- c(DNAoutput, paste(K, collapse=""))
    } else {
      DNAoutput <- rbind(DNAoutput, K)
    }
  }
  if (!is.null(filename)) {
    if (format=="nexus") {
      if (is.null(sampletime)) {
        genomenames <- paste("G", 1:length(sampleID), "_", sampleID, sep="")
      } else {
        genomenames <- paste("G", 1:length(sampleID), "_", sampletime, sep="")
      }
      DNAoutput <- cbind(genomenames, DNAoutput)
      write(t(DNAoutput), file=filename, sep="\t", append=TRUE, ncolumns=2)
      write(";\nend;", file=filename, append=TRUE)
    } else if (format=="fasta") {
      if (is.null(sampletime)) {
        write(paste(">lcl|G1_", sampleID[1], sep=""), file=filename)
      } else {
        write(paste(">lcl|G1_", sampletime[1], sep=""), file=filename)
      }
      k <- 1
      while (k<length(ref.strain)) {
        write(paste(DNAoutput[1,k:min(length(ref.strain),(k+49))], collapse=""), file=filename, append=TRUE)
        k <- k+50
      }
      if (ngenomes>1) {
        for (i in 2:ngenomes) {
          if (is.null(sampletime)) {
            write(paste(">lcl|G", i, "_", sampleID[i], sep=""), file=filename, append=TRUE)
          } else {
            write(paste(">lcl|G", i, "_", sampletime[i], sep=""), file=filename, append=TRUE)
          }
          k <- 1
          while (k<length(ref.strain)) {
            write(paste(DNAoutput[i,k:min(length(ref.strain),(k+49))], collapse=""), file=filename, append=TRUE)
            k <- k+50
          }
        }
      }
    }
  } else {
    return(DNAoutput)
  }
}
