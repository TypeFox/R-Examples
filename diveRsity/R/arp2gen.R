################################################################################
# arp2gen: arlequin file conversion to genepop
################################################################################
#' @export
arp2gen <- function(infile){
  # test if the file exists
  flForm <- strsplit(infile, split = "\\.")[[1]]
  if(substr(infile, 1, 2) == "./"){
    flForm <- flForm[-1]
  } else if(substr(infile, 1, 3) == "../"){
    flForm <- flForm[-(1:2)]
  }
  if(length(flForm) > 3){
    stop("There were multiple '.' characters in your file name!")
  }
  tstfile <- paste(flForm[1], ".gen", sep = "")
  # define a fastscan function
  if(!file.exists(tstfile)){
    fastScan <- function(fname) {
      s <- file.info(fname)$size
      buf <- readChar(fname, s, useBytes = TRUE)
      # replace Mac encoded line endings
      if(length(grep("\r", buf)) != 0L){
        buf <- gsub("\r", "\n", buf)
        buf <- gsub("\n\n", "\n", buf)
      }
      return(strsplit(buf, "\n", fixed = TRUE, useBytes = TRUE)[[1]])
    }
    
    # scan infile
    dat <- fastScan(infile)
    
    # strip needless whitespace
    dat <- gsub("^\\s+|\\s+$", "", dat)
    
    # some safeguards
    dataType <- grep("*datatype=*", tolower(dat))
    if(strsplit(dat[dataType], "=")[[1]][2] != "MICROSAT"){
      stop("Data are not in 'MICROSAT' format!")
    }
    
    # extract the relavant information (nloci, npops etc.)
    
    # missing data character
    missDataLine <- grep("*missingdata=*", tolower(dat))
    missData <- noquote(substr(dat[missDataLine],
                               nchar(dat[missDataLine]) - 1,
                               nchar(dat[missDataLine]) - 1))
    
    # samples sizes
    sampSizeLine <- grep("*samplesize=*", tolower(dat))
    if(length(sampSizeLine) > 1){
      sampNpos <- sapply(sampSizeLine, function(i){
        return(regexpr("=", dat[i])[1])
      })
    }
    popSizes <- as.numeric(substr(dat[sampSizeLine],
                                  start = sampNpos+1,
                                  stop = nchar(dat[sampSizeLine])))
    
    # number of population samples
    npops <- length(popSizes)
    
    # number of loci
    sampStrt <- grep("*sampledata=*", tolower(dat))
    
    # adjust sample starts for possible white space
    strts <- sapply(sampStrt, function(x){
      if(dat[(x+1)] == ""){
        return(x + 2)
      } else {
        return(x + 1)
      }
    })
    
    # define pop ends
    ends <- strts + ((popSizes * 2) - 1)
    
    nloci <- length(strsplit(dat[strts[1]], split = "\\s+")[[1]]) - 2
    
    # extract genotypes
    popGeno <- lapply(seq_along(strts), function(i){
      return(dat[strts[i]:ends[i]])
    })
    
    # check that popsizes are consistent
    popSzcheck <- sapply(popGeno, function(x) length(x)/2)
    if(!all(identical(popSzcheck, popSizes))){
      stop("Failed! Please make sure that your file is formatted correctly.")
    }
    
    # create a vector of odd indexes for each pop
    popIdx <- lapply(popGeno, function(x){
      return(seq(1, length(x), 2))
    })
    
    # paste alleles together
    popList <- lapply(seq_along(popGeno), function(i){
      al1 <- matrix(unlist(strsplit(popGeno[[i]][popIdx[[i]]], 
                                    split = "\\s+")), nrow = popSizes[i],
                    byrow = TRUE)[,-(1:2)]
      al2 <- matrix(unlist(strsplit(popGeno[[i]][(popIdx[[i]] + 1)], 
                                    split = "\\s+")), nrow = popSizes[i],
                    byrow = TRUE)
      tst <- matrix(paste(al1, al2, sep = ""), nrow = popSizes[i])
      tst <- cbind(paste(rep("pop", nrow(tst)), i, " ,", sep = ""), tst)
      # tidy up
      rm(al1, al2)
      z <- gc()
      rm(z)
      # replace missing data with genepop format
      if(nchar(tst[1,2]) == 4){
        tst[tst == paste(missData, missData, sep = "")] <- "0000"
      } else {
        tst[tst == paste(missData, missData, sep = "")] <- "000000"
      }
      out <- apply(tst, 1, function(x){
        return(paste(x, collapse = "\t"))
      })
      out <- c("POP", out)
      #     out <- rep(NA, nrow(tst))
      #     for(j in 1:nrow(tst)){
      #       out[j] <- paste(tst[j,], collapse = "\t")
      #     }
      # tidy up
      rm(tst)
      z <- gc()
      rm(z)
      return(out)
    })
    
    # A genepop file can not be written easily
    
    # Generate the outfile name
    outfile <- strsplit(infile, "\\.")[[1]]
    if(length(outfile) >= 2){
      outfile <- paste(outfile[-length(outfile)], collapse = ".")
    } else {
      outfile <- outfile[1]
    }
    
    # construct the file
    loci <- paste("locus", 1:nloci, sep = "")
    loci <- c(paste(outfile, "_gen_converted", sep = ""), loci)
    
    # outfile object
    of <- c(loci, unlist(popList))
    
    # define a file connection
    out <- file(paste(outfile, ".gen", sep = ""), "w")
    for(i in 1:length(of)){
      cat(of[i], "\n", file = out, sep = "")
    }
    close(out)
    return(TRUE)
  } else {
    return(NULL)
  }
}
################################################################################
# END - arp2gen
################################################################################