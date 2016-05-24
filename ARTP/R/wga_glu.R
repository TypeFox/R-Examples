# History Jul 02 2008 Initial coding
#         Nov 12 2008 Add glu.mergeAlleles
#         Dec 12 2008 Add glu.qcsummary
#         Dec 16 2008 Add options list to glu.mergeAlleles
#                     Add function glu.renameAlleles
#         Feb 02 2009 Update glu.qcsummary for samples
#                     Add function getMAF.control
#         Feb 17 2009 Update getMAF.control for subsetting the phenotype data
#         May 21 2009 Add glu.partition function
#         Jun 10 2009 Update glu.partition for cases and controls
#         Jul 07 2009 Add function glu.getCounts
#         Jul 20 2009 Fix bug in glu.partition with obtaining ccLevels
#         Jul 22 2009 Add error checks in glu.partition
#         Aug 10 2009 Add glu.r2 function
#         Oct 19 2009 Add glu.mergePhenoGeno function
#         Nov 02 2009 NO CHANGE
#         Dec 22 2009 Add glu.setUp function for 
#         Dec 28 2009 Put glu.setUp in wga_util2.R
#                     Add options for glu.partition
#         Jan 11 2010 Update temporary files in glu.partition 
#         Jan 14 2010 In glu.partition, allow data to be passed in
#         Sep 13 2010 Add function for obtaining r2 matrix
#         Oct 18 2010 Add function for getting the snps in genes
#         Oct 19 2010 Add function for retuning the number of bins using GLU
#         Oct 22 2010 Add function for returning the snps in bins based on
#                     an R^2 threshold and MAF.
#         Oct 25 2010 Change option names in glu.nBins
#         Nov 12 2010 Add function to convert to a ped file for Haploview
#         Nov 16 2010 Change options in ldmatrix and add function to add missing snps
#         Jan 21 2011 Convert mat to a matrix in addMissingSNPs

# Function to open a connection and transfrom data
glu.transform <- function(infile, inFormat=NULL, outFormat="ldat") {

  # infile      Path to file to open a connection
  #             No default
  # inFormat    GLU format of infile  
  #             The default is NULL
  # outFormat   The GLU format of the transformed data
  #             The default is "ldat".

  command <- "glu transform"
  if (!is.null(inFormat)) command <- paste(command, " -f ", inFormat, sep="")
  command <- paste(command, " ", infile, sep="")
  if (!is.null(outFormat)) command <- paste(command, " -F ", outFormat, sep="")
  
  fid <- pipe(command, open="r")

} # END: glu.transform

# Function for changing the alleles using GLU
glu.mergeAlleles <- function(snp1.list, snp2.list, temp.list=NULL,
                     op=NULL) {

  # snp1.list and snp2.list can also have the field "qcsummary".

  # snp1.list     Base genotype data list
  # snp2.list     
  # temp.list
  ###########################################################
  # op            List with names:
  # outfile       Updated genotype file for snp2.list
  # include.loci
  # flip2         0 or 1 for flipping alleles in snp2.list
  ###########################################################

  # Check the options list
  op <- default.list(op,
         c("outfile", "flip2"), list("out.txt", 0))

  outfile      <- op$outfile
  flip2        <- op$flip2
  include.loci <- getListName(op, "include.loci")
  rm(op)
  temp <- gc(verbose=FALSE)

  temp.list <- check.temp.list(temp.list)
  snp1.list <- check.snp.list(snp1.list)
  snp2.list <- check.snp.list(snp2.list)

  # Check if the qc.summary files exist
  tmp1 <- getListName(snp1.list, "qcsummary")
  if (is.null(tmp1)) {
    tmp1 <- getTempfile(temp.list$dir, prefix="snp1", ext=".txt")  
    temp <- glu.qcsummary(snp1.list, out.loci=tmp1, include.loci=include.loci)
  }
  tmp2 <- getListName(snp2.list, "qcsummary")
  if (is.null(tmp2)) {
    tmp2 <- getTempfile(temp.list$dir, prefix="snp2", ext=".txt") 
    temp <- glu.qcsummary(snp2.list, out.loci=tmp2, include.loci=include.loci)
  }

  # Get the snp names and alleles
  vars  <- c("LOCUS", "NUM_ALLELES", "ALLELES")
  tlist <- list(file=tmp1, file.type=3, delimiter="\t", header=1) 
  s1    <- getColumns(tlist, vars, temp.list=temp.list)
  tlist$file <- tmp2
  s2    <- getColumns(tlist, vars, temp.list=temp.list)
  names(s1$ALLELES)     <- s1$LOCUS
  names(s1$NUM_ALLELES) <- s1$LOCUS
  names(s2$ALLELES)     <- s2$LOCUS
  names(s2$NUM_ALLELES) <- s2$LOCUS

  snps <- intersect(s1$LOCUS, s2$LOCUS)
  
  # Remove snps = "*" (last row of qc.summary file)
  snps <- snps[snps != "*"]
  if (!length(snps)) {
    print("No intersecting SNPs")
    return(0)
  }
  a1   <- s1$ALLELES[snps]
  a2   <- s2$ALLELES[snps]
  n1   <- as.numeric(s1$NUM_ALLELES[snps])
  n2   <- as.numeric(s2$NUM_ALLELES[snps])

  # Delete temporary files
  if (temp.list$delete) {
    file.remove(tmp1)
    file.remove(tmp2)
  }
  
  rm(s1, s2, tmp1, tmp2)
  temp <- gc()

  # Check for more than 2 alleles
  temp <- (n1 > 2) + (n2 > 2)
  if (any(temp)) {
    stop("NUM_ALLELES > 2 for some SNPs")
  }

  # Seperate the alleles
  nsnp     <- length(a1)
  mat1     <- matrix(data="", nrow=nsnp, ncol=2)
  mat2     <- matrix(data="", nrow=nsnp, ncol=2)
  mat1[,1] <- substr(a1, 1, 1)
  mat1[,2] <- substr(a1, 3, 3)
  mat2[,1] <- substr(a2, 1, 1)
  mat2[,2] <- substr(a2, 3, 3)
  rm(a1, a2)
  temp <- gc()

  # Create a temporary file for the allele re-naming
  tmp <- getTempfile(temp.list$dir, prefix="allele", ext=".txt")  

  # Open and initialize the file
  temp <- c("SNP", "OLD_ALLELES", "NEW_ALLELES")
  fid  <- writeVecToFile(temp, tmp, colnames=NULL, type=3, close=0,
                           sep="\t")
  
  reverse <- 2:1

  # Loop over each snp
  for (i in 1:nsnp) {
    n2i <- n2[i]
    n1i <- n1[i]
    if (!n2i) next

    m1 <- mat1[i, 1:n1i]
    m2 <- mat2[i, 1:n2i]

    # See if the alleles match
    flag <- (all(m2 %in% m1)) || (all(m1 %in% m2))

    # If they agree and there is no flip, then do nothing
    if (flag) {
      if (!flip2) next 

      # Alleles agree and there is a flip
      if (n2i == 2) {
        old <- paste(m2, collapse=",", sep="")
        new <- paste(m2[reverse], collapse=",", sep="")
      } else {
        old <- m2[1]
        if (n1i == 2) {
          # Get the other allele
          if (old == m1[1]) {
            new <- m1[2]
          } else {
            new <- m1[1]
          }
        } else if (n1i == 1) {
          new <- m1[1]
          if (new == old) next
        } else {
          # n1[i] is 0. Do nothing
          next
        }
      }
    } else {
      if (!n1i) next

      # Alleles do not match
      old <- paste(m2, collapse=",", sep="")
      new <- paste(m1, collapse=",", sep="")
      if (n1i == 1) {
        new <- m1[1]
        old <- m2[1]
      } else if (n2i == 1) {
        new <- m1[1]
      }
      if (old == new) next 
    }

    # Write to file
    temp <- c(snps[i], old, new)
    write(temp, file=fid, ncolumns=3, sep="\t")

  } 

  close(fid)

  # Rename the alleles
  temp <- glu.renameAlleles(snp2.list, tmp, outfile=outfile)

  # Delete temporary file
  if (temp.list$delete) file.remove(tmp)
    
  0

} # END: glu.mergeAlleles

# Function for calling ginfo
glu.ginfo <- function(snp.list, outloci="out.txt") {

  snp.list  <- check.snp.list(snp.list)
  type      <- snp.list$file.type
  if (type == 2) type <- "ldat"
  if (type == 3) type <- "sdat"
  
  snp     <- paste(snp.list$dir, snp.list$file, sep="")
  command <- paste("glu ginfo --outputloci=", outloci, sep="")
  command <- paste(command, " -f ", type, " ", snp, ":unique=n", sep="")
  ret     <- callOS(command)
  ret

} # END: glu.ginfo

# Function to get the common snps
glu.commonSNPs <- function(slist, temp.list=NULL, outfile=NULL) {

  # slist     List of sublists of type snp.list

  temp.list <- check.temp.list(temp.list)
  tmpfile   <- getTempfile(temp.list$dir, prefix="snp", ext=".txt") 
  tlist     <- list(file=tmpfile, file.type=3, delimiter="\t", header=1) 
  n         <- length(slist)
  
  for (i in 1:n) {
    temp <- glu.ginfo(slist[[i]], outloci=tmpfile)
    temp <- getColumns(tlist, "LOCUS", temp.list=temp.list)
    if (i == 1) {
      snps <- temp[[1]]
    } else {
      snps <- intersect(snps, temp[[1]])
    }
  }
  snps <- unique(snps)

  if (!is.null(outfile)) write(snps, file=outfile, ncolumns=1)

  snps

} # END: glu.commonSNPs

# Function to call qc.summary
glu.qcsummary <- function(snp.list, out.loci="out.txt", include.loci=NULL,
                  include.samples=NULL) {

  snp.list  <- check.snp.list(snp.list)
  
  snp     <- paste(snp.list$dir, snp.list$file, sep="")
  command <- paste("glu qc.summary --locusout=", out.loci, sep="")
  command <- paste(command, " ", snp, sep="")
  if (!is.null(include.loci)) {
    command <- paste(command, " --includeloci=", include.loci, sep="")
  }
  if (!is.null(include.samples)) {
    command <- paste(command, " --includesamples=", include.samples, sep="")
  }
  ret <- callOS(command)
  ret

} # END: glu.qcsummary

# Function to rename alleles
glu.renameAlleles <- function(snp.list, alleleFile, outfile="out.txt") {

  snp.list <- check.snp.list(snp.list)
  type     <- snp.list$file.type
  if (type == 2) type <- "ldat"
  if (type == 3) type <- "sdat"
  
  snp     <- paste(snp.list$dir, snp.list$file, sep="")
  command <- paste("glu transform -o ", outfile, sep="")
  command <- paste(command, " -f ", type, sep="")
  command <- paste(command, " --renamealleles=", alleleFile, sep="")
  command <- paste(command, " ", snp, sep="")
  ret     <- callOS(command)
  ret


} # END: glu.renameAlleles

# Function to compute the MAF in controls using GLU
getMAF.control <- function(snp.list, pheno.list, data, op=NULL, temp.list=NULL) {
 
  # Check the lists
  snp.list   <- check.snp.list(snp.list)
  pheno.list <- check.pheno.list(pheno.list)
  temp.list  <- check.temp.list(temp.list)

  # Get the subset of controls
  x <- read.table(pheno.list$file, header=pheno.list$header, sep=pheno.list$delimiter)

  temp <- getListName(pheno.list, "subsetData")
  if (!is.null(temp)) {
    x <- subsetData.list(x, temp)
  }

  vars <- c(pheno.list$id.var, pheno.list$response.var)
  x    <- removeOrKeepCols(x, vars, which=1)
  x    <- unfactor.all(x)

  # Get the controls 
  temp <- x[, pheno.list$response.var] == 0
  x    <- removeOrKeepRows(x, temp, which=1) 
  
  # Create a temporary file for the subject ids
  tmp <- getTempfile(temp.list$dir, prefix="samples", ext=".txt") 
  
  # Write the samples ids
  write(x[, pheno.list$id.var], file=tmp, ncolumns=1)

  # Create a temporary file for the snp ids
  tmp2 <- getTempfile(temp.list$dir, prefix="snps", ext=".txt") 
  
  # Write the samples ids
  write(data[, "SNP"], file=tmp2, ncolumns=1)

  rm(x, temp, vars)
  temp <- gc(verbose=FALSE)

  # Create a temporary file for GLU output
  gluout <- getTempfile(temp.list$dir, prefix="glu", ext=".txt") 

  # Call GLU
  temp <- glu.qcsummary(snp.list, out.loci=gluout, include.loci=tmp2,
                  include.samples=tmp)

  # Add the MAF column
  temp <- list(file=gluout, file.type=3, delimiter="\t", header=1, id.var="LOCUS", vars="MAF", names="MAF")
  data <- addColumn(data, "SNP", temp)

  # Delete files
  if (temp.list$delete) {
    file.remove(tmp)
    file.remove(gluout)
    file.remove(tmp2)
  }

  data

} # END: getMAF.control

# Function to get info on which snps belong to which study
glu.partition <- function(snp.list, pheno.list, temp.list=NULL, op=NULL) {

  # Check the lists
  op         <- default.list(op, c("byCC"), list(0), error=c(0))
  snp.list   <- check.snp.list(snp.list)
  pheno.list <- check.pheno.list(pheno.list)
  temp.list  <- check.temp.list(temp.list)
  pheno.list <- default.list(pheno.list, c("partition.var", "is.the.data"),
                  list("ERROR", 0), error=c(1, 0))
  idVal      <- temp.list$id
  
  if (temp.list$dir %in% c("", " ")) temp.list$dir <- "./"

  # Check if byCC option was specified
  byCC <- op$byCC
  if (byCC) {
    temp <- pheno.list[["cc.var", exact=TRUE]]
    if (is.null(temp)) {
      stop("ERROR: cc.var not specified in pheno.list")
    } 
  } else {
    ccLevels <- 1
  }
 
  pvar <- pheno.list$partition.var
  id   <- pheno.list$id.var
  x    <- readTable(pheno.list$file, pheno.list)
  pv   <- makeVector(unfactor(x[, pvar]))
  id   <- makeVector(unfactor(x[, id]))
  if (byCC) {
    cv <- makeVector(unfactor(x[, pheno.list$cc.var]))
    ccLevels <- unique(cv)
  }
  nr <- nrow(x)
  rm(x)
  gc()
 
  # Get the partition levels
  levels <- unique(pv)
  
  # Create a temporary file for the subject ids
  tmpfile <- getTempfile(temp.list$dir, prefix=paste("samples_", idVal, sep=""), ext=".txt") 

  # Create a temporary file for the qc output
  tmpqc <- getTempfile(temp.list$dir, prefix=paste("qcout_", idVal, sep=""), ext=".txt") 

  # GLU variables
  vars <- c("NUM_GENOTYPES", "GENOTYPE_COUNTS")

  # Loop over the levels
  x <- NULL
  for (lev in levels) {
    print(lev)
    temp0 <- (pv == lev)
    temp0[is.na(temp0)] <- FALSE
    
    for (cc in ccLevels) {

      if (byCC) {
        temp <- temp0 & (cv == cc)
        temp[is.na(temp)] <- FALSE
      } else {
        temp <- temp0
      }
      s <- sum(temp)
      print(s)
      if (!s) next
      
      # Write out the ids
      write(id[temp], file=tmpfile, ncolumns=1)

      # Call qc.summary
      temp <- glu.qcsummary(snp.list, out.loci=tmpqc, include.loci=NULL,
                  include.samples=tmpfile)

      # New variable names
      vname <- paste(pvar, ".", lev, sep="")
      new  <- c(vname, paste(vname, ".COUNTS", sep=""))
      if (byCC) new <- paste(new, ".", cc, sep="")

      # Add the info
      if (is.null(x)) {
        tlist <- list(what="character", returnMatrix=1, include.row1=1, delimiter="\t")
        x <- scanFile(tmpqc, tlist)
        x <- x[, c("LOCUS", "NUM_GENOTYPES", "GENOTYPE_COUNTS")]
        x <- data.frame(x)
        x <- unfactor.all(x)
        x <- renameVar(x, "LOCUS", "SNP")
        x <- renameVar(x, vars[1], new[1])
        x <- renameVar(x, vars[2], new[2])
        x[, new[1]] <- as.numeric(x[, new[1]])
        partVars <- new[1]
      } else {
        tlist <- list(file=tmpqc, file.type=3, header=1, delimiter="\t",
                      id.var="LOCUS", vars=vars, names=new, type=c("N", "C"))
        x <- addColumn(x, "SNP", tlist)
        partVars <- c(partVars, new[1])
      }

      # Delete files
      file.remove(tmpfile)
      file.remove(tmpqc)

    } # END: for (cc in ccLevels) 
  } # END: for (lev in levels

  # Determine which snps belonged to all studies
  v <- "ALL_PARTS"
  x[, v] <- 1
  for (var in partVars) x[, v] <- as.numeric(x[, v] & (x[, var] > 0)) 

  out <- op[["outfile", exact=TRUE]]
  if (!is.null(out)) writeTable(x, out)

  # Return the SNPs
  temp <- x[, v] == 1
  temp[is.na(temp)] <- FALSE
  snps <- makeVector(x[temp, "SNP"])

  snps

} # END: glu.partition 

# Function to get (study specific geno counts)
glu.getCounts <- function(snp.list, pheno.list, temp.list=NULL, op=NULL) {

  # op        List with names
  #  byCC     0 or 1
  #           The default is 1
  #  by.var   Name of variable on phenotype data for seperate counts
  #           The default is NULL
  #  snps     List of type file.list or a character vector of snps
  #           The default is NULL

  # Check lists
  snp.list   <- check.snp.list(snp.list)
  pheno.list <- check.pheno.list(pheno.list)
  temp.list  <- check.temp.list(temp.list)
  op         <- default.list(op, c("byCC"), list(1))
  
  # Get the by variable
  byVar <- getListName(op, "by.var")
  if (is.null(byVar)) {
    byFlag <- 0
  } else {
    byFlag <- 1
  }
  
  # Get the case-control variable
  ccVar <- getListName(pheno.list, "cc.var")
  if (is.null(ccVar)) op$byCC <- 0
  ccFlag <- op$byCC

  # Get the snps
  snpFlag <- 1
  snps <- getListName(op, "snps")
  if (is.null(snps)) snpFlag <- 0
  if (snpFlag) {
    if (typeof(snps) != "character") {
      snps <- scan(snps$file, what="character", sep="\n")
    }

    # Write these snps out to a temporary file
    tmpfile <- getTempfile(temp.list$dir, prefix="snps", ext=".txt") 
    write(snps, file=tmpfile, ncolumns=1)
    rm(snps)
    gc()
  }
  
  # Create temporary file for GLU output
  tmpglu <- getTempfile(temp.list$dir, prefix="gluOut", ext=".txt") 

  # Create string 
  str <- paste("glu transform -F ldat -o ", tmpglu, sep="")
  if (snpFlag) str <- paste(str, " --includeloci=", tmpfile, sep="")
  str <- paste(str, " ", snp.list$file, sep="")
  print(str)
  temp <- callOS(str)

  # Update snp.list
  snp.list$file        <- tmpglu
  snp.list$file.type   <- 2
  snp.list$stream      <- 0
  pheno.list$keep.vars <- c(pheno.list$id.var, ccVar, byVar)

  # Get the data
  tlist <- list(include.row1=0, include.snps=0, return.type=1, MAF=0,
                missing=1, snpNames=1, orderByPheno=1, return.pheno=1)
  temp  <- try(getData.1(snp.list, pheno.list, temp.list, op=tlist),
               silent=TRUE)
  if (class(temp) == "try-error") {
    print(temp)
    stop("ERROR loading data")
  }

  snpData   <- temp$data
  #missing   <- temp$missing
  snpNames  <- temp$snpNames
  delimiter <- "\t"
  nsnps     <- length(snpData)

  # Get the phenotype data
  phenoData.list <- temp$phenoData.list
  phenoData0     <- phenoData.list$data
  nr             <- nrow(phenoData0)

  if (ccFlag) cc <- makeVector(phenoData0[, ccVar])
  if (byFlag) by <- makeVector(phenoData0[, byVar])
  rm(temp, phenoData.list, phenoData0)
  gc()

  # Remove temporary files
  file.remove(tmpfile)
  file.remove(tmpglu)

  # Get the unique levels of ccVar and byVar
  if (ccFlag) {
    ccLevels <- unique(cc)
    ncc      <- length(ccLevels)
  } else {
    ncc      <- 1
  }
  if (byFlag) {
    byLevels <- unique(by)
    nby      <- length(byLevels)
  } else {
    nby      <- 1
  }

  # Initialize return data frame
  ret <- data.frame(snpNames)
  colnames(ret) <- "SNP"
  rm(snpNames)
  gc()

  # Loop over each level
  for (i in 1:nby) {
    if (byFlag) {
      print(byLevels[i])
      tempBy <- by == byLevels[i]
      tempBy[is.na(tempBy)] <- FALSE
    } else {
      tempBy <- rep(TRUE, times=nr)
    }
    for (j in 1:ncc) {
      if (byFlag) {
        print(ccLevels[j])
        tempCC <- cc == ccLevels[j]
        tempCC[is.na(tempCC)] <- FALSE
      } else {
        tempCC <- rep(TRUE, times=nr)
      }

      # Get the correct subset
      rows <- tempBy & tempCC
      
      # Get the new variable name
      if (byFlag) new <- byLevels[i]
      if (ccFlag) {
        new <- paste(new, ".", ccLevels[j], sep="")
      }
      ret[, new] <- NA 

      # Loop over each SNP
      for (k in 1:nsnps) {
        snp         <- as.integer(getVecFromStr(snpData[k], delimiter=delimiter))
        counts      <- getGenoCounts(snp[rows], exclude=c(NA, NaN), check=1) 
        ret[k, new] <- paste(counts, collapse="|", sep="")
      }
    }
  } # END: for (i in 1:nby)

  # Output file
  out <- getListName(op, "outfile")
  if (!is.null(out)) writeTable(ret, out)

  ret

} # END: glu.getCounts

# Function to get r2 and dprime.
# Subjects to use is determined by pheno.list$subsetData. Set pheno.list
#  to use all subjects in the genotype data.
# Use snp.list$snpNames to specify the snps
glu.r2 <- function(snp.list, pheno.list=NULL, op=NULL) {

  # op            List with names
  #  outfile

  # Check lists
  snp.list   <- check.snp.list(snp.list)
  if (!is.null(pheno.list)) {
    pheno.list <- check.pheno.list(pheno.list)
    pflag      <- 1
  } else {
    pflag      <- 0
  }
  op         <- default.list(op, c("outfile"), list("ERROR"), error=1)
  temp.list  <- getListName(op, "temp.list")
  temp.list  <- check.temp.list(temp.list)
  
  # Write out subject ids
  if (pflag) {
    temp  <- getPhenoData(pheno.list, temp.list=temp.list)$data
    tfile <- getTempfile(temp.list$dir, prefix="subs_r2_", ext=".txt")
    temp  <- makeVector(temp[, pheno.list$id.var])
    write(temp, file=tfile, ncolumns=1)
  }
 
  # Write out snps
  snps <- getListName(snp.list, "snpNames")
  sflag <- !is.null(snps)
  if (sflag) {
    snpfile <- getTempfile(temp.list$dir, prefix="SNP_r2_", ext=".txt")
    write(snp.list$snpNames, file=snpfile, ncolumns=1)
    snpStr <- paste(" --includeloci=", snpfile, sep="")
  } else {
    snpStr <- paste(" --includeloci=", snp.list$snpNames.list$file, sep="")
  }
  
  # Call GLU
  str <- "glu tagzilla --skipbinning -r 0"
  str <- paste(str, " --saveldpairs=", op$outfile, sep="")
  str <- paste(str, snpStr, sep="")
  if (pflag) str <- paste(str, " --includesamples=", tfile, sep="")
  str <- paste(str, " ", snp.list$file, sep="")
  callOS(str)
  
  # Delete file
  if (temp.list$delete) {
    if (sflag) file.remove(snpfile) 
    if (pflag) file.remove(tfile)
  }

  # Read in and edit the glu output
  x <- read.table(op$outfile, header=1, sep="\t", as.is=FALSE)
  x <- unfactor.all(x)
  temp <- x[, "LNAME1"] == x[, "LNAME2"]
  temp[is.na(temp)] <- FALSE
  x <- removeOrKeepRows(x, temp, which=-1)

  writeTable(x, op$outfile)

  x

} # END: glu.r2

# Function to merge phenotype and genotype data
glu.mergePhenoGeno <- function(snp.list, pheno.list, op=NULL) {

  # op          List of options for mergePhenoGeno and temp.list
  #  outfile
  #  which      0-2 for the coding

  snp.list   <- check.snp.list(snp.list)
  pheno.list <- check.pheno.list(pheno.list)
  temp.list  <- getListName(op, "temp.list")
  temp.list  <- check.temp.list(temp.list)
  snp.list   <- default.list(snp.list, "snpNames", list("ERROR"), error=1)

  # Temp file for ldat and snpNames
  tmpfile <- getTempfile(temp.list$dir, prefix="SNP", ext=".ldat")
  tmp     <- getTempfile(temp.list$dir, prefix="SNPNAMES", ext=".txt")
  temp    <- getListName(snp.list, "snpNames")
  write(temp, file=tmp, ncolumns=1)  

  str <- "glu transform -F ldat -o "
  str <- paste(str, tmpfile, sep="")
  str <- paste(str, " --includeloci=", tmp, sep="")
  str <- paste(str, " --orderloci=", tmp, sep="")
  str <- paste(str, " ", snp.list$file, sep="") 
  callOS(str)

  snp.list$file      <- tmpfile
  snp.list$file.type <- 2
  snp.list$delimiter <- "\t"
  x <- mergePhenoGeno(snp.list, pheno.list, temp.list=temp.list, op=op)

  if (temp.list$delete) {
    file.remove(tmpfile)
    file.remove(tmp)
  }

  x

} # END: glu.mergePhenoGeno

# Function to add missing SNPs in ld matrix
addMissingSNPs <- function(snps, mat) {

  n <- length(snps)
  nr <- nrow(mat)
  if (n == nr) return(mat)

  cnames <- colnames(mat)
  temp   <- !(snps %in% cnames)
  miss   <- snps[temp]
  nmiss  <- length(miss) 
  mat    <- as.matrix(mat)
  colnames(mat) <- NULL
  rownames(mat) <- NULL
  mat    <- rbind(mat, matrix(0, nrow=nmiss, ncol=nr))
  mat    <- cbind(mat, matrix(0, nrow=nrow(mat), ncol=nmiss))
  cnames <- c(cnames, miss)  
  rownames(mat) <- cnames
  colnames(mat) <- cnames
  diag(mat)     <- 1
  
  mat
}

# Function for returning an ld matrix
glu.ldMatrix <- function(snps, genoFile, op=NULL) {

  temp.list <- op[["temp.list", exact=TRUE]]
  temp.list <- check.temp.list(temp.list)
  op <- default.list(op, c("maxdist", "minmaf", "measure", "order", "addMiss"), 
                 list(200, 0.05, "r2", 1, 1))

  # Determine where it is running
  #host <- callOS("hostname", inter=TRUE)
  #if (host == "biowulf.nih.gov) {
  #  bioFlag <- 1
  #} else {
  #  bioFlag <- 0
  #}

  snps <- unique(snps)
  nsnps <- length(snps)

  if (nsnps == 1) {
    x <- matrix(data=1.0, nrow=1, ncol=1)
    rownames(x) <- snps
    colnames(x) <- snps
    return(x)
  }

  # Write snps to a temporary file
  tempfile <- getTempfile(temp.list$dir, prefix=paste("_tmplist_", temp.list$id, "_", sep=""), ext=".txt")
  write(snps, file=tempfile, ncolumns=1)

  # Temporary output file
  tempout <- getTempfile(temp.list$dir, prefix=paste("_tmpout_", temp.list$id, "_", sep=""), ext=".txt")
  
  # Create the glu call command
  str <- paste("glu ld.matrix --maxdist=", op$maxdist, " --minmaf=", op$minmaf, 
               " --measure=", op$measure, " --includeloci=", tempfile,
               " --output=", tempout, sep="")
  # Loci description file
  #temp <- op[["loci.file", exact=TRUE]]
  #if ((is.null(temp)) && ())
 
  str <- paste(str, " ", genoFile, sep="")
  callOS(str)
 
  # Read in the matrix
  x <- read.table(tempout, as.is=TRUE, header=1, row.names=1, sep="\t")
  n <- ncol(x)
  for (i in 1:(n-1)) {
    cols <- (i+1):n
    x[i, cols] <- x[cols, i]
  }
  temp <- is.na(x)
  x[temp] <- 0

  # Delete temp files
  if (temp.list$delete) {
    file.remove(tempfile)
    file.remove(tempout)
  }

  # Add missing snps if needed
  if (op$addMiss) x <- addMissingSNPs(snps, x)

  # Check for error
  if (nsnps != nrow(x)) stop("ERROR in glu.ldMatrix: incorrect dimension of LD matrix")

  # Order the same as the input list of snps
  if (op$order) x <- x[snps, snps]
  
  x

} # glu.ldMatrix

# Function for getting the snps in a list of genes
glu.snps_in_genes <- function(gene.list, op=NULL) {

  keep <- c("LOCUS",  "CHROMOSOME",  "LOCATION", "FEATURE_NAME", "FEATURE_TYPE")

  op <- default.list(op, c("d", "u", "keep.vars"), list(70000, 70000, keep))
  temp.list <- op[["temp.list", exact=TRUE]]
  temp.list <- check.temp.list(temp.list)
  dir <- temp.list$dir  
  id <- temp.list$id

  genoFile <- op[["geno.file", exact=TRUE]]
  genoFlag <- !is.null(genoFile)

  tmpfile <- getTempfile(dir, paste("genedb_", id, sep=""), ext=".txt")
  tmpgene <- getTempfile(dir, paste("genes_", id, sep=""), ext=".txt")
  write(gene.list, file=tmpgene, ncolumns=1)

  str <- paste("glu genedb.find_snps -u ", op$u, " -d ", op$d, 
               " -o ", tmpfile, " ", tmpgene, sep="")
  print(str)
  callOS(str)

  x <- loadData.table(tmpfile)
  keep <- op[["keep.vars", exact=TRUE]]
  if (!is.null(keep)) x <- x[, keep]

  if (genoFlag) {
    str <- paste("glu ginfo --outputloci=", tmpfile, ":c=1 ", genoFile, sep="")
    print(str)
    callOS(str)
    snps <- scan(tmpfile, what="character", sep="\n")
    temp <- x[, "LOCUS"] %in% snps
    x <- x[temp, ]
    rm(snps)
    gc()
  }

  if (temp.list$delete) {
    file.remove(tmpfile)
    file.remove(tmpgene)
  }

  temp <- op[["out.file", exact=TRUE]]
  if (!is.null(temp)) writeTable(x, temp)
  
  x

} # END: glu.snps_in_genes

# Function for returning the number of bins
glu.nBins <- function(snps, genoFile, op=NULL) {

  op <- default.list(op, c("r2.threshold", "min.MAF", "temp.list", "return.data", "snpsIsFile", "samplesIsFile"), 
                    list(0.8, 0.05, list(), 0, 0, 0))
  temp.list <- check.temp.list(op$temp.list)

  snpFlag <- op$snpsIsFile
  if (!snpFlag) {
    snps <- unique(snps)
    nsnps <- length(snps)
    if (nsnps == 1) return(1)

    # Write snps to a temporary file
    tempfile <- getTempfile(temp.list$dir, prefix=paste("_tmplist_", temp.list$id, "_", sep=""), ext=".txt")
    write(snps, file=tempfile, ncolumns=1)
  } else {
    tempfile <- snps
  }

  samples  <- op[["samples", exact=TRUE]]
  sflag    <- !is.null(samples)
  sampFlag <- op$samplesIsFile
  if (sflag) {
    if (!sampFlag) {
      # Write samples to a temporary file
      tempsamp <- getTempfile(temp.list$dir, prefix=paste("_tmpsamp_", temp.list$id, "_", sep=""), ext=".txt")
      write(samples, file=tempsamp, ncolumns=1)
    } else {
      tempsamp <- samples
    }
  }

  # Temporary output file
  tempout <- getTempfile(temp.list$dir, prefix=paste("_tmpout_", temp.list$id, "_", sep=""), ext=".txt")
  
  # Create the glu call command
  str <- paste("glu ld.tagzilla -r ", op$r2.threshold, " -a ", op$min.MAF, 
               " --includeloci=", tempfile,
               " -O ", tempout, sep="")
  if (sflag) str <- paste(str, " --includesamples=", tempsamp, sep="")
  str <- paste(str, " ", genoFile, sep="")
  callOS(str)
 
  # Read in the matrix
  x <- read.table(tempout, as.is=TRUE, header=1, sep="\t")
  n <- max(as.numeric(x[, "BINNUM"]), na.rm=TRUE)
  
  # Delete temp files
  if (temp.list$delete) {
    if (!snpFlag) file.remove(tempfile)
    file.remove(tempout)
    if ((sflag) && (!sampFlag)) file.remove(tempsamp)
  }

  if (op$return.data) return(x)

  n

} # END: glu.nBins

# Function for returning the snps in a set of bins
glu.SNP_bins <- function(snps, genoFile, op=NULL) {

  # op    
  #  samples        Character vector of sample ids

  op <- default.list(op, c("return.data"), list(1))

  x <- glu.nBins(snps, genoFile, op=op) 

  # LNAME         LOCATION    MAF           BINNUM      DISPOSITION

  # For each bin, choose the snp with the largest MAF
  bin <- as.numeric(makeVector(x[, "BINNUM"]))
  snp <- makeVector(x[, "LNAME"])  
  maf <- as.numeric(makeVector(x[, "MAF"]))
  rm(x)
  gc()

  ubins <- unique(bin)
  nbins <- length(ubins)
  ret   <- character(nbins)
  for (i in 1:nbins) {
    temp <- bin == ubins[i]
    if (sum(temp) == 1) {
      ret[i] <- snp[temp]
    } else {
      j <- which.max(maf[temp])
      ret[i] <- (snp[temp])[j]   
    }
  }

  ret

} # END: glu.SNP_bins

# Function to convert to a ped file
glu.create_ped <- function(snp.list, pheno.list, op=NULL) {
 
  # snp.list        File must be a format GLU can read
  # pheno.list      With names response.var, gender.var, male, female
  # op       
  #   map.type      0=file from GLU, 1=format for haploview

  snp.list <- check.snp.list(snp.list)
  pheno.list <- check.pheno.list(pheno.list)
  pheno.list <- default.list(pheno.list, 
                c("response.var", "gender.var", "male", "female"), 
                list("ERROR", "ERROR", "MALE", "FEMALE"), 
                error=c(1, 1, 0, 0))
  op <- default.list(op, c("temp.list", "map.type"), list(list(), 0))
  temp.list <- check.temp.list(op$temp.list)

  x         <- loadData.table(pheno.list)
  vars      <- c(pheno.list$id.var, pheno.list$response.var, pheno.list$gender.var)
  x         <- removeOrKeepCols(x, vars)
  subids    <- makeVector(x[, pheno.list$id.var])
  gender0   <- makeVector(x[, pheno.list$gender.var])
  response0 <- makeVector(x[, pheno.list$response.var])
  gender    <- integer(length(gender0))
  temp      <- gender0 %in% pheno.list$male
  gender[temp] <- 1
  temp      <- gender0 %in% pheno.list$female
  gender[temp] <- 2
  response  <- integer(length(response0))
  temp     <- response0 == 0
  temp[is.na(temp)] <- FALSE
  response[temp] <- 1
  temp     <- response0 == 1
  temp[is.na(temp)] <- FALSE
  response[temp] <- 2

  rm(x, gender0, response0)
  gc()

  # See if snps are in a file already
  snp.file <- op[["snp.file", exact=TRUE]]
  snpFlag  <- !is.null(snp.file) 

  # Create file for sample ids
  dir <- temp.list$dir
  id  <- temp.list$id
  tmpids <- getTempfile(dir, prefix=paste("ids_", id, "_", sep=""), ext=".txt")
  write(subids, file=tmpids, ncolumns=1)

  # Temporary output file
  tmpfile <- getTempfile(dir, prefix=paste("out_", id, "_", sep=""), ext=".ped")

  str <- paste("glu transform -F ped -o ", tmpfile, sep="")
  str <- paste(str, " --includeloci=", snp.file, sep="")
  str <- paste(str, " --includesamples=", tmpids, sep="")
  str <- paste(str, " ", snp.list$file, sep="")
  print(str)
  callOS(str)

  # File has columns: Family id, subject id, father id, mother id, gender, phenotype, genotypes
  ped <- scan(tmpfile, what="character", sep="\n")
  n   <- length(ped)
  
  for (i in 1:n) {
    x      <- getVecFromStr(ped[i], delimiter=" ")
    row    <- match(x[2], subids)
    x[5:6] <- c(gender[row], response[row]) 
    ped[i] <- paste(x, collapse=" ", sep="")  
  }

  out <- op[["out.ped", exact=TRUE]]
  if (!is.null(out)) write(ped, file=out, ncolumns=1)
  rm(ped)
  gc() 

  out <- op[["out.map", exact=TRUE]]
  if (!is.null(out)) {
    # Get the map file name
    len <- nchar(tmpfile)
    tmpmap <- substr(tmpfile, 1, len-4)
    tmpmap <- paste(tmpmap, ".map", sep="")
    file.copy(tmpmap, out, overwrite=TRUE)
  }

  if (op$map.type) {
    x <- scan(out, what="character", sep=" ")
    x <- matrix(x, byrow=TRUE, ncol=4)
    x <- x[, c(2, 4)]
    write.table(x, file=out, sep=" ", row.names=FALSE, quote=FALSE, col.names=FALSE)
  }

  if (temp.list$delete) {
    file.remove(tmpids)
    file.remove(tmpfile)
    file.remove(tmpmap)
  }

}

