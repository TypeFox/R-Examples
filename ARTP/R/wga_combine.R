# History: Mar 10 2008  Initial coding
#          Apr 08 2008  Add combine.snp to combine biowulf results
#                       and compute permutation p-values for each snp
#          Jun 06 2008  Update locusMap.list options
#          Jun 25 2008  Call getMap; do not factor the chromsome vector
#                       Put combine.rhs.tests options in a list
#          Jul 03 2008  Add getSNPlists
#          Sep 19 2008  Add option to appendFiles
#          Oct 21 2008  Fix bug in appendFiles
#          Oct 27 2008  Bug fixes for combining studies functions
#          Nov 03 2008  Add meta.analysis function
#          Nov 14 2008  Allow snps to be NULL in meta analysis
#          Nov 14 2008  Add getResults.scan function
#          Nov 21 2008  Allow header to be output in appendFiles 
#          Dec 10 2008  Remove getResults.scan
#                       Generalize the major.allele options in
#                       meta.analysis.
#          Dec 18 2008  Add variables to data frame in meta.analysis
#          Jan 09 2009  Add globalAlleles function
#          Jan 13 2009  Add copy.vars option to meta.analysis
#          Jan 16 2009  Add MAF to output in global.alleles function
#          Jan 21 2009  Generalize meta.analysis function to read in
#                       the global minor alleles.
#          Jan 30 2009  Add allele.flip function
#          Feb 10 2009  Use only controls in global.allele function
#          Feb 12 2009  In meta.analysis, add option to compute method 2.
#                       Add meta.mima.scan function for type 3 meta-analysis
#          Mar 05 2009  Move unused code to wga_unused.R
#          Apr 01 2009  Check for errors in LOR.flip
#          Apr 03 2009  Fix bug in LOR.flip (take intersection of snps)
#          May 26 2009  Add option to use variances instead of se in meta-analysis.
#          Jun 03 2009  Allow snp names to be passed into rows option
#                       in meta.mima.scan function.
#          Jun 04 2009  In meta.mima.scan add option to pass in mods
#          Jun 11 2009  Update getNperm function for permuations in 
#                       seperate files
#          Jul 08 2009  Make meta.analysis, allele.flip functions more efficient
#                       reading in the data. 
#          Aug 17 2009  Use stringsAsFactors option in allele.flip, lor.flip, meta.analysis
#          Aug 18 2009  Change in LOR.flip to only include snps passed in allele vector
#          Aug 19 2009  globalAlleles: add doc, remove options, check for errors
#                       Use stringsAsFactors = FALSE option
#          Aug 20 2009  Change in appendFiles for file.type=3 and vars=NULL
#          Nov 10 2009  Return number of rows in file in appendFiles()
#          Jan 14 2010  Add option to appendFiles to check for the correct number
#                       of lines in each file read in.
#          Oct 14 2010  Combine meta anlaysis functions into function meta
#          Dec 09 2010  Change output varible names in meta.anlaysis
#                       Add OR and CI
#          Jan 12 2011  Remove fun.op in combObserved.study
#                       Add the mima_v2 function
#          Jan 31 2011  Add function to convert obs file to data frame
#                       Remove getListName calls

# Function to combine the times
# Returns a list of total time and average time in minutes
combineTime <- function(dir="./", pattern="GenABEL_time", outfile=NULL) {

  print(pattern)

  # Get the files
  files <- list.files(path=dir, pattern=pattern)

  sum <- 0
  n   <- 0

  for (file in files) {
    t00 <- NULL
    load(paste(dir, file, sep=""))

    sum <- sum + t00
    n   <- n + 1 
  }

  # Get the time in minutes
  sum <- sum/60

  print("Total time:")
  print(sum)
  print("Average time:")
  print(sum/n)

  t00 <- list(total=sum, avg=sum/n)

  if (!is.null(outfile)) save(t00, file=outfile)

  t00

} # END: getTime

# Function to append files 
appendFiles <- function(dir, op=NULL) {

  # op              List
  #   nlines.list   List of type file.list with file names and number of lines
  #                 in each file
  #                 The default is NULL

  # Assumes common file.type, header, delimiter among all files
  op <- default.list(op, 
         c("outfile", "header", "file.type", "delimiter", "print.files"),
         list("ERROR", 1, 3, "\t", 1), error=c(1, 0, 0, 0, 0))  

  error <- 0
  files <- op[["file", exact=TRUE]]
  if (is.null(files)) files <- list.files(path=dir, pattern=op$pattern)
  vars    <- op[["vars", exact=TRUE]]
  varFlag <- !is.null(vars)
  missingFiles <- NULL
  
  # For option nlines
  nlist <- op[["nlines.list", exact=TRUE]]
  nFlag <- 0
  if (!is.null(nlist)) {
    nlist <- check.file.list(nlist)
    nlist <- default.list(nlist, c("id.var", "nlines.var"), 
                         list("ERROR", "ERROR"), error=c(1, 1))
    temp  <- loadData(nlist$file, nlist)
    nl.f  <- as.character(temp[, nlist$id.var])
    nl.n  <- as.integer(temp[, nlist$nlines.var])
    names(nl.n) <- nl.f
    nFlag <- 1

    # See if any files are missing
    temp <- !(nl.f %in% files)
    if (any(temp)) {
      missingFiles <- nl.f[temp]
      print(missingFiles)
      print("WARNING in appendFiles: Above files were not found")
      error <- 1
    }

    rm(temp, nlist, nl.f)
    gc()
  }

  # Open the output file
  fid <- file(op$outfile, open="w")

  row1     <- op[["row1", exact=TRUE]]
  row1Flag <- !is.null(row1)
  op$row1  <- NULL
  sep0     <- op$delimiter
    
  i         <- 1
  dir       <- checkForSep(dir)
  nlines    <- 0
  op$stream <- 0
  if (varFlag) {
    op$file.type    <- 3
    op$method       <- 2
    op$returnMatrix <- 1
    op$what         <- "character"
    op$include.row1 <- 1
    ncolumns        <- length(vars)
    delimiter       <- op$delimiter

    # Get the number of columns in a data set
    temp    <- paste(dir, files[1], sep="")
    op$ncol <- getNcols(temp, list(delimiter=delimiter)) 
    
    if (!row1Flag) write(vars, file=fid, ncolumns=ncolumns, sep=delimiter)
  } else {
    header       <- op$header 
    op$delimiter <- "\n"
    if (op$file.type == 3) op$file.type=2
  }    

  if (row1Flag) {
    write(row1, file=fid, ncolumns=length(row1), sep=sep0)
    rm(row1, sep0, row1Flag)
    temp <- gc()
  }

  errorFiles  <- NULL
  print.files <- op$print.files
  for (file in files) {
    if (print.files) print(file)
    temp <- paste(dir, file, sep="")
    temp <- loadData(temp, op)
    if (varFlag) {
      # temp is a matrix
      if (nFlag) {
        if (nl.n[file] != nrow(temp)) errorFiles <- c(errorFiles, file)
      }

      temp <- temp[, vars]
      write(t(temp), file=fid, ncolumns=ncolumns, sep=delimiter)
      nlines <- nlines + nrow(temp)
    } else {
      if (nFlag) {
        if (nl.n[file] != length(temp)) errorFiles <- c(errorFiles, file)
      }

      if ((i > 1) && (header)) temp <- temp[-1]
      write(temp, file=fid, ncolumns=1)
      nlines <- nlines + length(temp)
      i <- i + 1
    }
  }
  close(fid)

  if (!is.null(errorFiles)) {
    print(errorFiles)
    print("WARNING in appendFiles: Above files have the incorrect number of rows")
    error <- 1
  }

  temp <- paste("nlines in output file = ", nlines, sep="")
  print(temp)

  list(error=error, errorFiles=errorFiles, nlines=nlines, missingFiles=missingFiles)

} # END: appendFiles

# Function to return a list of the observed vectors in each study
getObserved.study <- function(study.list) {

  obs <- list()

  for (i in 1:length(study.list)) {
    temp     <- study.list[[i]]$obs.list
    obs[[i]] <- getObserved.vec(temp) 
  } 
  obs

} # END: getObserved.study

# Function to return the observed vector
getObserved.vec <- function(obs.list) {

  # Set up the list of options
  tlist <- list(returnMatrix=1, start.row=1, method=2,
               include.row1=1, what=double(0), stop.row=-1) 

  tlist$delimiter <- obs.list$delimiter
  tlist$file.type <- obs.list$file.type

  # Load the data
  temp <- loadData(obs.list$file, tlist)    
  temp <- makeVector(temp)
 
  temp

} # END: getObserved.vec

# Function to get the SNP names in each study
getSNPnames.study <- function(obs.list) {

  snames <- list()
  for (i in 1:length(obs.list)) {
    snames[[i]] <- names(obs.list[[i]])
  }
  snames

} # END: getSNPnames.study

# Function to get the snps for each gene and the intersection for each gene
geneSNPs.study <- function(snames, snps, genes, method=2) {

  # snames      List of character vectors
  # snps        SNPs from gene file
  # ugenes

  snps.list <- list()
  nstudy    <- length(snames)
  for (i in 1:nstudy) snps.list[[i]] <- list()
  csnps.list <- list()
  ugenes <- unique(genes)
  
  for (gg in ugenes) { 
    # Get the snps for each gene
    temp  <- (genes == gg)
    snp2  <- snps[temp]

    # Initialize the intersection
    csnps.list[[gg]] <- snp2

    for (i in 1:nstudy) {
      temp <- snp2 %in% snames[[i]]
      temp <- snp2[temp]
      snps.list[[i]][[gg]] <- temp

      # Get the intersection of the snps
      csnps.list[[gg]] <- intersect(csnps.list[[gg]], temp)       
    }
  } 

  # Change the snps for each gene, depending on the method
  if (method == 2) {
    for (i in 1:nstudy) snps.list[[i]] <- csnps.list
  }

  # Create vector of unique snps names for a combined study
  csnps <- NULL
  for (i in 1:length(csnps.list)) {
    csnps <- c(csnps, csnps.list[[i]])
  }

  list(snps.list=snps.list, csnps.list=csnps.list, csnps=csnps)

} # END: geneSNPs.study

# Function to get the number of rows to read in
getReadN.double <- function(nc, maxMB=100) {

  return(maxMB*1000*1000/(8*nc))

} # END: getReadN.double

# Function to check input study lists
check.study.list <- function(study.list) {

  for (i in 1:length(study.list)) {
    temp <- study.list[[i]]
    temp <- default.list(temp,
           c("obs.list", "perm.list"),
           list("ERROR", "ERROR"), error=c(1, 1))

    # Check the observed list
    temp$obs.list <- check.obs.list(temp$obs.list)

    # Check the permutation list
    temp$perm.list <- check.obs.list(temp$perm.list)

    # Check the data list
    dlist <- temp[["data.list", exact=TRUE]]
    if (!is.null(dlist)) {
      dlist <- check.data.list(dlist)
      temp$data.list <- dlist
    }

    study.list[[i]] <- temp
  }

  study.list

} # END: check.study.list

# Function to check the data list
check.data.list <- function(data.list) {

  data.list <- default.list(data.list, 
      c("file", "file.type", "delimiter", "header", "snp.var", "vars"),
      list("ERROR", 3, " ", 1, "SNP", "ERROR"),
      error=c(1, 0, 0, 0, 0, 1))

  data.list
  
} # END: check.data.list

# Function to check the observed list
check.obs.list <- function(obs.list) {

  obs.list <- default.list(obs.list, 
      c("file", "file.type", "delimiter"),
      list("ERROR", 3, " "),
      error=c(1, 0, 0))

  obs.list
  
} # END: check.obs.list

# Function to check the permutation list
check.perm.list <- function(perm.list) {

  perm.list <- default.list(perm.list, 
      c("file", "file.type", "delimiter"),
      list("ERROR", 3, " "),
      error=c(1, 0, 0))

  perm.list
  
} # END: check.perm.list

# Function to check the gene list
check.gene.list <- function(gene.list) {

  gene.list <- default.list(gene.list,
   c("file", "delimiter", "file.type", "gene.var", "snp.var", "header"),
   list("ERROR", "\t", 3, "Gene", "SNP", 1),
   error=c(1, 0, 0, 0, 0, 0))

  gene.list
  
} # END: check.gene.list

# Function to check the group list
check.group.list <- function(group.list) {

  group.list <- default.list(group.list,
   c("file", "delimiter", "file.type", "group.var", "gene.var", "header"),
   list("ERROR", "\t", 3, "Category", "Gene", 1),
   error=c(1, 0, 0, 0, 0, 0))

  group.list
  
} # END: check.group.list

# Function to get the combined observed values
combObserved.study <- function(obs, csnps, fun, op=NULL) {

  # obs     List of observed vectors

  comb.obs <- NULL
  for (i in 1:length(obs)) {
   
    # Get the snps we need
    temp <- obs[[i]][csnps]

    comb.obs <- fun(comb.obs, temp, op=op)
  } 
  
  # Call the combine function again
  comb.obs <- fun(comb.obs, NULL, op=op)
  names(comb.obs) <- csnps

  comb.obs 

} # END: combObserved.study

# Function to get the combined permutation values
combPerm.study <- function(fid.list, tlist, csnps, cnames, nperm,
                  total.read, fun, fun.op=NULL) {

  comb.perm <- NULL
 
  for (i in 1:length(fid.list)) {
    # Read the data
    temp <- readPermData(fid.list[[i]], tlist[[i]], csnps, cnames[[i]],
                         nperm=nperm, total.read=total.read[i])

    if (!length(temp)) return(NULL)
    total.read[i] <- temp$total.read

    # Combine 
    comb.perm <- fun(comb.perm, temp$data, op=fun.op)
  }  
 
  # Call the combine function again
  comb.perm <- fun(comb.perm, NULL, op=fun.op)
  colnames(comb.perm) <- csnps

  list(comb.perm=comb.perm, total.read=total.read)

} # END: combPerm.study

# Function to read permutation data
readPermData <- function(fid, tlist, csnps, cnames, nperm=Inf,
                  total.read=0) {

  temp <- scanFile(fid, tlist)

  if (!length(temp)) return(NULL)
  colnames(temp) <- cnames 
  total.read <- total.read + nrow(temp)
  
  temp <- subsetPermData(temp, csnps, total.read=total.read, 
                           nperm=nperm)

  list(data=temp, total.read=total.read)

} # END: readPermData

# Function to get the correct subset of the perm data
subsetPermData <- function(data, csnps, total.read=0, nperm=Inf) {

  # Get the snps we need
  data <- removeOrKeepCols(data, csnps, which=1)

  if (total.read >= nperm) {
    remove <- total.read - nperm

    nr <- nrow(data)
       
    # Remove the last n rows
    if (remove) data <- removeOrKeepRows(data, 1:(nr-remove), which=1)
  }

  data 

} # END: subsetPermData

# Function to compute the number of snps for each study
getNcols.study <- function(study.list) {

  for (i in 1:length(study.list)) {
    perm.list <- study.list[[i]]$perm.list
    nc        <- perm.list[["nsnps", exact=TRUE]]
    if (is.null(nc)) {
      perm.list$nsnps <- getNcols(perm.list$file, perm.list)
      study.list[[i]]$perm.list <- perm.list
    }
  }

  study.list

} # END: getNcols.study

# Function to compute the number of permutations for each study
getNperm.study <- function(study.list) {

  for (i in 1:length(study.list)) {
    perm.list <- study.list[[i]]$perm.list
    perm.list$nperm <- getNperm(perm.list)
    study.list[[i]]$perm.list <- perm.list
  }

  study.list

} # END: getNrows.study

# Function to get the number of permutations
getNperm <- function(perm.list) {

  nr <- perm.list[["nperm", exact=TRUE]]
  if (!is.null(nr)) return(nr)

  # See if the directory option was specified
  dir <- perm.list[["dir", exact=TRUE]]
  if (!is.null(dir)) {
    pattern <- perm.list[["pattern", exact=TRUE]]
    files   <- list.files(path=dir, pattern=pattern)
    dir     <- checkForSep(dir)
    files   <- paste(dir, files, sep="")
  } else {
    files <- perm.list$file
  }
  type <- perm.list$file.type

  nr <- 0
  for (f in files) {
    temp <- getNrows(f, file.type=type)
    nr   <- nr + temp - 1
  }

  nr

} # END: getNperm

# Get the minimum number of permutations among the studies
getMinPerm.study <- function(study.list) {

  nperm <- Inf
  for (i in 1:length(study.list)) {
    perm.list <- study.list[[i]]$perm.list
    temp      <- perm.list[["nperm", exact=TRUE]]
    nperm     <- min(temp, nperm)
  }
  nperm

} # END: getMinPerm.study

# Function to compute gene-level test stats for a vector
geneStats.init <- function(vec, ugenes, csnps.list, fun=which.max) {

  # vec             Vector with snps as names
  # ugenes
  # csnps.list      List of vectors of snp names for each gene in ugenes
  # fun             which.max or which.min  (for now)

  n <- length(ugenes)

  # Set up a data frame
  temp <- c("gene", "n.snp", "most.sig.snp", "test.stat")
  ret  <- initDataFrame(n, c("c", "n", "c", "n"), rownames=ugenes,
                        colnames=temp)

  for (gg in ugenes) {
    snp2                    <- csnps.list[[gg]]
    tests                   <- vec[snp2]
    temp                    <- fun(tests)
    test                    <- tests[temp]
    sigsnp                  <- snp2[temp] 
    if (!length(sigsnp)) {
      sigsnp <- ""
      test   <- NA
    }
    ret[gg, "gene"]         <- gg
    ret[gg, "n.snp"]        <- length(snp2)
    ret[gg, "most.sig.snp"] <- sigsnp
    ret[gg, "test.stat"]    <- test
  } 

  ret

} # END: geneStats.init

# Function to compute gene level test stat for observed results
geneStat.obs <- function(vec, ugenes, csnps.list, fun=NULL, ...) {

  ret        <- rep(NA, times=length(vec))
  names(ret) <- ugenes

  for (gg in ugenes) {
    snp2    <- csnps.list[[gg]]
    tests   <- vec[snp2]
    ret[gg] <- fun(tests, ...)
  } 

  ret

} # END: geneStat.obs

# Function to compute gene level test stat for permutation results
geneStat.perm <- function(perm.mat, obs.vec, count.vec, ugenes,
                       csnps.list, ret.mat=0, pvalues=1) {
  
  # perm.mat    Matrix of p-values
  # obs.vec     Vector of observed p-values

  if (ret.mat) {
    temp.mat <- matrix(data=NA, nrow=nrow(perm.mat),
                       ncol=length(ugenes))
    colnames(temp.mat) <- ugenes
  } else {
    temp.mat <- NULL
  }

  if (pvalues) {
    fun <- min
  } else {
    fun <- max
  }

  for (gg in ugenes) {
    snp2 <- csnps.list[[gg]]
 
    # Get the subset of the matrix
    temp <- removeOrKeepCols(perm.mat, snp2, which=1)

    # Get the test stat for each permutation
    temp <- apply(temp, 1, fun, na.rm=TRUE)

    # Add up the counts
    if (pvalues) {
      count.vec[gg] <- count.vec[gg] +
                      sum(temp <= obs.vec[gg], na.rm=TRUE)
    } else {
      count.vec[gg] <- count.vec[gg] +
                      sum(temp >= obs.vec[gg], na.rm=TRUE)
    }

    if (ret.mat) temp.mat[, gg] <- temp 
  } 

  list(count.vec=count.vec, temp.mat=temp.mat)

} # END: geneStat.perm

# Function to return a list of fids
getFID.study <- function(study.list, name="perm.list") {

  fid <- list()
  for (i in 1:length(study.list)) {
    temp     <- getListName(study.list[[i]], name)
    file     <- temp[["file", exact=TRUE]]
    type     <- temp[["file.type", exact=TRUE]]
    temp     <- list(file.type=type, open="r")
    fid[[i]] <- getFID(file, temp)
  }
  fid

} # END: getFID.study

# Function to combine pvalues
combPvalue.log <- function(combined, newData, op=NULL) {

  if (!is.null(newData)) {
    temp <- -2*log(newData)
  } else {
    temp <- 0
  }
  if (!is.null(combined)) temp <- temp + combined
  return(temp)

} # END: combPvalue.log

# Function to return additional variables needed to combine values
combVars.study <- function(study.list, csnps, temp.list=NULL) {

  ret <- list()
  for (i in 1:length(study.list)) {
    dlist <- getListName(study.list[[i]], "data.list")
    if (!is.null(dlist)) {
      vars <- c(dlist$snp.var, dlist$vars)
      temp <- getColumns(dlist, vars, temp.list=temp.list)
      for (var in dlist$vars) {
        names(temp[[var]]) <- temp[[dlist$snp.var]]
      }
      for (var in vars) temp[[var]] <- temp[[var]][csnps]
      for (var in dlist$vars) temp[[var]] <- as.numeric(temp[[var]])
      ret[[i]] <- temp
    }
  }
  ret
  
} # END: combVars.study

# Function to compute a meta analysis between 2 or more studies.
# If some SNP results were based on the major allele and some on the 
#  minor allele, then set major.allele, and the combined results
#  will be based on the minor allele.
meta.analysis <- function(study.list, op=NULL) {

  # study.list         List of sublists with the folowing fields
  #  file              Name of file containing results.
  #                    No default.
  #  file.type         1, 3, or 4.
  #                    The default is 3.
  #  header            0 or 1
  #                    The default is 1.
  #  delimiter         The default is "\t"
  #  snp.var           Name of the variable in file for the SNPs.
  #                    The default is "SNP".
  #  lor.var           NULL or name of the log(odds) variable
  #                    The default is NULL.
  #  lor.se.var        NULL or name of the variable for the standard error
  #                    of log(odds).
  #                    The default is NULL.
  #  lor.ucl.var       NULL or the name of the variable for the upper confidence
  #                    limit for log(odds).
  #                    The default is NULL.
  #  lor.conf          Confidence level for lor.ucl.var.
  #                    The default is 0.95. 
  #  or.var
  #  or.se.var
  #  or.ucl.var
  #  or.conf
  #  flip.lor          0 if no flipping the lor is to be done, 1 if all snps are
  #                    to be flipped, or a
  #                    list containing the fields "var", "operator",
  #                    and "value" to specify which snps are based to
  #                    be flipped.
  #                    Example: list(var="MAF", operator=">", value=0.85) specifies
  #                    that the SNPs with MAF > 0.85 are to be flipped.
  #                    The default is 0.
  #  copy.vars         Variables to copy onto the output data set.
  #                    These variables names are not prefixed with the
  #                    study name.
  #                    The default is NULL.
  #  copy2.vars        Variables to copy onto the output data set.
  #                    These variables names are prefixed with the
  #                    study name.
  #                    The default is NULL.
  #  allele.var        Name of the allele variable to compare to the global
  #                    allele (see allele.list in op).
  #                    If the allele does not match the global allele, then
  #                    the log-odds ratio wil be flipped.
  #                    The default is NULL.
  #  n.case            (For method=2)
  #  n.control         (For method=2)
  ################################################################
  # op            List with the following fields
  #  snps         Character vector of snps for the analysis
  #               If NULL, then all common snps are used
  #               The default is NULL
  #  outfile      Output file.
  #               The default is NULL
  #  add.vars     Variables to appear in outfile for each study.
  #               ex. c("lor.var", "lor.se.var", "or.var")
  #               The default is NULL
  #  add.names    Variable names for add.vars.
  #               ex. c("LOR", "LOR.SE", "OR")
  #               The default is NULL
  #  method       1 or 2
  #               The default is 1.
  ########################################################################
  #  allele.list  List of type file.list that contains the info about
  #               the global alleles. This list takes precedence over flip.lor.
  #    allele.var Name of the allele variable to compare for each study.
  #               No default.
  #    snp.var    SNP variable
  #               No default.
  ########################################################################

  nstudy <- length(study.list)

  op <- default.list(op, c("method"), list(1), error=c(0))
  method <- op$method

  # Check for global allele list
  temp <- op[["allele.list", exact=TRUE]]
  if (!is.null(temp)) {
    temp <- default.list(temp, c("allele.var", "snp.var"), 
            list("ERROR", "ERROR"), error=c(1, 1)) 
    alleleFlag <- 1
  } else {
    alleleFlag <- 0
  }

  # Check each list
  for (i in 1:nstudy) {
    tlist <- default.list(study.list[[i]],
     c("file", "file.type", "header", "delimiter", "snp.var", 
       "lor.conf", "or.conf", "flip.lor", "name", "allele.var"),
     list("ERROR", 3, 1, "\t", "SNP", 0.95, 0.95, 0, 
          paste("study", i, sep=""), "TEMP"), 
     error=c(1, 0, 0, 0, 0, 0, 0, 0, 0, alleleFlag))

    lor.flag     <- !is.null(tlist[["lor.var", exact=TRUE]])
    lor.se.flag  <- !is.null(tlist[["lor.se.var", exact=TRUE]])
    lor.ucl.flag <- !is.null(tlist[["lor.ucl.var", exact=TRUE]])
    or.flag      <- !is.null(tlist[["or.var", exact=TRUE]])
    or.se.flag   <- !is.null(tlist[["or.se.var", exact=TRUE]])
    or.ucl.flag  <- !is.null(tlist[["or.ucl.var", exact=TRUE]])

    if (!lor.flag && !or.flag) {
      stop("One of lor.var or or.var must be specified")
    }
    if (!lor.flag) {
      # OR se must be given
      if (!or.se.flag && !or.ucl.flag) {
        stop("One of or.se.var or or.ucl.var must be specified")
      }
      tlist$lor.var <- "LOR_TEMPVAR1234"
    } 

    if (lor.se.flag) {
      lor.ucl.flag      <- 0
      tlist$lor.ucl.var <- NULL
      tlist$lor.conf    <- NULL
    } else {
      tlist$lor.se.var  <- "LORSE_TEMPVAR1234"
    }

    if (!or.flag) {
      # LOR se must be given
      if (!lor.se.flag && !lor.ucl.flag) {
        stop("One of lor.se.var or lor.ucl.var must be specified")
      }
      tlist$or.var <- "OR_TEMPVAR1234"
    } 

    if (or.se.flag) {
      or.ucl.flag      <- 0
      tlist$or.ucl.var <- NULL
      tlist$or.conf    <- NULL
    } else {
      tlist$or.se.var  <- "ORSE_TEMPVAR1234"
    }

    tlist$test.var   <- "TEST_TEMPVAR1234"
    tlist$pvalue.var <- "PVALUE_TEMPVAR1234"
    tlist$flip.var   <- "FLIP_TEMPVAR1234"

    if (!lor.ucl.flag) tlist$lor.conf <- NULL
    if (!or.ucl.flag)  tlist$or.conf  <- NULL

    # Check for flip var
    temp <- tlist[["flip.lor", exact=TRUE]]
    if (is.list(temp)) {
      temp <- default.list(temp, c("var", "operator", "value"),
               list("ERROR", "ERROR", "ERROR"), error=c(1, 1, 1))
      mvar <- temp$var
    } else {
      mvar <- NULL
    }


    # Check for the number of cases and controls if method = 2
    if (method == 2) {
      if (is.null(tlist[["n.case", exact=TRUE]]) || is.null(tlist[["n.case", exact=TRUE]])) {
        stop("ERROR: n.case/n.control not specified")
      } 
    }

    # Get the variables to keep
    temp <- mvar
    if (lor.flag)    temp  <- c(temp, tlist$lor.var)
    if (lor.se.flag) temp  <- c(temp, tlist$lor.se.var)
    if (lor.ucl.flag) temp <- c(temp, tlist$lor.ucl.var)
    if (or.flag)      temp <- c(temp, tlist$or.var)
    if (or.se.flag)   temp <- c(temp, tlist$or.se.var)
    if (or.ucl.flag)  temp <- c(temp, tlist$or.ucl.var)
    tlist$numVars <- temp
    if (alleleFlag)   temp <- c(temp, tlist$allele.var)
    temp <- c(tlist$snp.var, temp)
    cvars <- tlist[["copy.vars", exact=TRUE]]
    if (!is.null(cvars)) temp <- c(temp, cvars)
    cvars <- tlist[["copy2.vars", exact=TRUE]]
    if (!is.null(cvars)) temp <- c(temp, cvars)
    rm(cvars)

    tlist$keepVars <- unique(temp)

    study.list[[i]] <- tlist
  } # END: for (i in 1:nstudy)

  snps <- op$snps

  # Read in the global alleles
  if (alleleFlag) {
    print("Reading global allele file")
    temp     <- op[["allele.list", exact=TRUE]]
    data     <- getColumns(temp, c(temp$snp.var, temp$allele.var))
    g.snp    <- data[[temp$snp.var]]
    g.allele <- data[[temp$allele.var]]
    names(g.allele) <- g.snp
    
    if (is.null(snps)) {
      snps <- g.snp
    } else {
      snps <- intersect(snps, g.snp)
    }
    rm(g.snp, data)
    temp <- gc()
  }

  # Read in the data, and then subset by snps and vars
  print("Reading in the data")
  for (i in 1:nstudy) {
    tlist <- study.list[[i]]
    print(tlist$name)
    if (tlist$file.type %in% c(3, 6, 8)) {
      tlist$method       <- 2
      tlist$what         <- "character"
      tlist$include.row1 <- tlist$header
    }
    temp  <- loadData(tlist$file, tlist)
    temp  <- removeOrKeepCols(temp, tlist$keepVars, which=1)
    if (!is.null(snps)) {
      ids   <- temp[, tlist$snp.var] %in% snps
      temp  <- removeOrKeepRows(temp, ids, which=1)
      rm(ids)
      gc()
    }
    tlist$data <- temp

    # Check for missing snps
    temp <- is.na(match(op$snps, temp[, tlist$snp.var]))
    if (any(temp)) {
      print(op$snps[temp])
      temp <- paste("The above SNPs were not found for study ", i, sep="")
      print(temp)
    }

    # Get the intersection of the snps for the analysis
    if (i == 1) {
      if (is.null(snps)) {
        snps <- tlist$data[, tlist$snp.var]
      } else {
        snps <- intersect(snps, tlist$data[, tlist$snp.var])
      }
    } else {
      snps <- intersect(snps, tlist$data[, tlist$snp.var])
    }

    study.list[[i]] <- tlist
  }
 
  # Order each data frame
  print("Ordering the data")
  for (i in 1:nstudy) {
    tlist <- study.list[[i]]
    temp  <- tlist$data
    rownames(temp) <- temp[, tlist$snp.var]
    temp  <- removeOrKeepRows(temp, snps, which=1)
    temp  <- data.frame(temp, stringsAsFactors=FALSE)
    temp  <- unfactor.all(temp)
    for (v in tlist$numVars) temp[, v] <- as.numeric(temp[, v])
    tlist$data <- temp
    study.list[[i]] <- tlist
  }

  # Check the alleles
  if (alleleFlag) {
    print("Checking the alleles")
    g.allele <- g.allele[snps]
    for (i in 1:nstudy) {
      tlist <- study.list[[i]]
      data  <- tlist$data
      data[, tlist$flip.var] <- 0
      temp <- g.allele != data[, tlist$allele.var]
      if (any(temp)) {
        data[temp, tlist$flip.var] <- 1
        tlist$flip.lor <- list(var=tlist$flip.var, operator="==", value=1)
      }
      tlist$data <- data
      study.list[[i]] <- tlist
    }
    rm(g.allele, data)
  }

  rm(tlist)
  temp <- gc()

  # Set up each data frame
  print("Setting up each data frame")
  for (i in 1:nstudy) {
    study.list[[i]] <- meta.setup(study.list[[i]])
  }

  # Do the meta analysis
  print("Starting meta-analysis")
  print(paste("Using method ", method, sep=""))
  if (method == 2) {
    # Compute the effective sample sizes
    totalN <- 0
    for (i in 1:nstudy) {
      tlist  <- study.list[[i]]
      temp   <- 2/(1/tlist$n.case + 1/tlist$n.control)
      totalN <- totalN + temp 
      study.list[[i]]$effN <- temp
    }
    # Compute the test statistic
    comb.test <- 0
    for (i in 1:nstudy) {
      tlist     <- study.list[[i]]
      data      <- tlist$data
      comb.test <- comb.test + sqrt(tlist$effN/totalN)*(data[, tlist$lor.var]/data[, tlist$lor.se.var])
    }
    comb.pval <- 2*pnorm(abs(comb.test), lower.tail=FALSE)
    comb <- data.frame(snps, comb.test, comb.pval, stringsAsFactors=FALSE)
    colnames(comb) <- c("SNP", "Meta.Test", "Meta.Pvalue")
    rm(comb.test, comb.pval, totalN, data)
    temp <- gc()
  } else if (method %in% c("1a", "1A")) {
    # Get the combined estimates
    comb.est <- 0
    comb.var <- 0
    for (i in 1:nstudy) {
      tlist    <- study.list[[i]]
      data     <- tlist$data
      # Compute variance
      data[, tlist$lor.se.var] <- data[, tlist$lor.se.var]*data[, tlist$lor.se.var]

      comb.est <- comb.est + data[, tlist$lor.var]/data[, tlist$lor.se.var]
      comb.var <- comb.var + 1/data[, tlist$lor.se.var]
    }
    comb.est  <- comb.est/comb.var
    comb.se   <- 1/sqrt(comb.var)
    comb.test <- comb.est/comb.se
    comb.pval <- 2*pnorm(abs(comb.test), lower.tail=FALSE)
    comb <- data.frame(snps, comb.est, comb.se, comb.test, comb.pval, stringsAsFactors=FALSE)
    colnames(comb) <- c("SNP", "Meta.Estimate", "Meta.SE", "Meta.Test", "Meta.Pvalue")
    rm(comb.est, comb.se, comb.test, comb.pval, comb.var)
    temp <- gc()
  } else {
    # Get the combined estimates
    comb.est <- 0
    comb.se  <- 0
    for (i in 1:nstudy) {
      tlist    <- study.list[[i]]
      data     <- tlist$data
      comb.est <- comb.est + data[, tlist$lor.var]/data[, tlist$lor.se.var]
      comb.se  <- comb.se + 1/data[, tlist$lor.se.var]
    }

    comb.est  <- comb.est/comb.se
    comb.se   <- sqrt(nstudy)/comb.se
    comb.test <- comb.est/comb.se
    comb.pval <- 2*pnorm(abs(comb.test), lower.tail=FALSE)
    comb <- data.frame(snps, comb.est, comb.se, comb.test, comb.pval, stringsAsFactors=FALSE)
    colnames(comb) <- c("SNP", "Meta.Estimate", "Meta.SE", "Meta.Test", "Meta.Pvalue")
    rm(comb.est, comb.se, comb.test, comb.pval)
    temp <- gc()
  }

  # Add OR and CI for method != 2
  if (method != 2) {
    beta <- as.numeric(comb[, "Meta.Estimate"])
    se   <- as.numeric(comb[, "Meta.SE"])
    comb[, "Meta.OR"] <- exp(beta)
    l    <- exp(beta - 1.96*se)
    l    <- round(l, digits=4)
    u    <- exp(beta + 1.96*se)
    u    <- round(u, digits=4)
    comb[, "Meta.OR.95CI"] <- paste("(", l, ", ", u, ")", sep="") 
    rm(beta, se, l, u)
    gc()
  }
  

  # Add copy.vars
  for (i in 1:nstudy) {
    tlist <- study.list[[i]]
    temp  <- tlist[["copy.vars", exact=TRUE]]
    if (!is.null(temp)) {
      temp <- removeOrKeepCols(tlist$data, temp, which=1)
      comb <- cbind(comb, temp)
    }
  }

  # Add copy2.vars
  for (i in 1:nstudy) {
    tlist <- study.list[[i]]
    old   <- tlist[["copy2.vars", exact=TRUE]]
    vars  <- paste(tlist$name, ".", old, sep="")
    if (!is.null(old)) {
      temp <- removeOrKeepCols(tlist$data, old, which=1)
      for (j in 1:length(old)) temp <- renameVar(temp, old[j], vars[j])  
      comb <- cbind(comb, temp)
    }
  }
  rm(old, vars)
  temp <- gc()

  # Add variables to the data frame for each study
  vars <- op[["add.vars", exact=TRUE]]
  add  <- op[["add.names", exact=TRUE]]
  if (!is.null(vars)) {
    if (length(vars) != length(add)) add <- vars
    for (i in 1:nstudy) {
      tlist <- study.list[[i]]
      data  <- tlist$data
      name  <- tlist$name
      j     <- 1
      for (var in vars) {
        temp <- paste(name, ".", add[j], sep="")
        comb[, temp] <- data[, tlist[[var]]]
        j <- j + 1
      }
    }
  } # END: if (!is.null(vars)) 

  # Sort by Pvalue
  comb <- sort2D(comb, "Meta.Pvalue")

  if (!is.null(op$outfile)) {
    write.table(comb, file=op$outfile, sep="\t", quote=FALSE, row.names=FALSE)
  }

  comb 

} # END: meta.analysis

# Function to change LOR, OR if based on a major allele
meta.update <- function(study) {

  temp <- study[["flip.lor", exact=TRUE]]
  if (is.null(temp)) return(study)
  flag <- is.list(temp) 
  if (!flag) {
    if (!temp) return(study)
  }

  # Get the rows to update
  if (flag) {
    rows <- subsetData.var(study$data, temp$var, temp$operator, 
                         temp$value, which=1, returnRows=1)
    if (!any(rows)) return(study)
  } else {
    rows <- 1:nrow(study$data)
  }
  data <- study$data

  # Check for LOR
  lorFlag <- study$lor.var %in% colnames(data)
  if (lorFlag) {
    data[, study$lor.var] <- as.numeric(data[, study$lor.var])
    lor <- data[, study$lor.var]
  }

  # Change LOR to -LOR, standard errors remain the same 
  if (lorFlag) {
    # See if ucl is given
    if (!is.null(study[["lor.ucl.var", exact=TRUE]])) {
      # Compute the standard error
      temp  <- study[["lor.conf", exact=TRUE]]
      temp  <- temp + (1-temp)/2
      zcrit <- qnorm(temp, lower.tail=TRUE)
      temp  <- data[, study$lor.ucl.var] - data[, study$lor.var]
      data[, study$lor.se.var] <- temp/zcrit
      study$lor.ucl.var <- NULL
      study$lor.conf    <- NULL
    }

    # Update the rows
    data[rows, study$lor.var] <- -data[rows, study$lor.var]
  }

  # Change OR to 1/OR, standard errors change 
  if (study$or.var %in% colnames(data)) {
    # See if ucl and lor.se are given
    flag1 <- is.null(study[["or.ucl.var", exact=TRUE]])
    flag2 <- study$lor.se.var %in% colnames(data)
    if (!flag1 && !flag2) {
      # Compute the standard error for LOR
      temp  <- study[["or.conf", exact=TRUE]]
      temp  <- temp + (1-temp)/2
      zcrit <- qnorm(temp, lower.tail=TRUE)
      if (lorFlag) {
        # Use original LOR or original OR
        temp <- log(data[, study$or.ucl.var]) - lor
      } else {
        temp <- log(data[, study$or.ucl.var]) - log(data[, study$or.var])
      }
      data[, study$lor.se.var] <- temp/zcrit
      study$or.ucl.var  <- NULL
      study$or.conf     <- NULL
      study$lor.ucl.var <- NULL
      study$lor.conf    <- NULL
    }

    # Change the standard error for 1/OR
    if (study$or.se.var %in% colnames(data)) {
      data[rows, study$or.se.var] <- data[rows, study$or.se.var]/data[rows, study$or.var]
    }

    # Now change OR
    data[rows, study$or.var] <- 1/data[rows, study$or.var]
  }

  study$data <- data
  study

} # END: meta.update

# Function to set up a data frame for a meta-analysis
meta.setup <- function(study) {

  study <- default.list(study,
           c("data", "lor.var", "lor.se.var", "or.var", "or.se.var",
             "lor.conf", "or.conf", "test.var", "pvalue.var"),
           list("ERROR", "LOR", "LOR.SE", "OR", "OR.SE", 0.95, 0.95,
                "TEST", "PVALUE"),
           error=c(1, 0, 0, 0, 0, 0, 0, 0, 0))
  temp <- study[["lor.ucl.var", exact=TRUE]]
  if (is.null(temp)) study$lor.conf <- NULL
  temp <- study[["or.ucl.var", exact=TRUE]]
  if (is.null(temp)) study$or.conf <- NULL

  # Call meta.update for major alleles
  study <- meta.update(study)

  data <- study$data
  study$data <- NULL
  gc()

  # Compute LOR  
  if (!(study$lor.var %in% colnames(data))) {
    data[, study$lor.var] <- log(data[, study$or.var])
  }

  # Compute OR  
  if (!(study$or.var %in% colnames(data))) {
    data[, study$or.var] <- exp(data[, study$lor.var])
  }

  # If confidence level is given compute standard errors
  temp <- study[["lor.conf", exact=TRUE]]
  if (!is.null(temp)) {
    temp  <- temp + (1-temp)/2
    zcrit <- qnorm(temp, lower.tail=TRUE)
    temp  <- data[, study$lor.ucl.var] - data[, study$lor.var]
    data[, study$lor.se.var] <- temp/zcrit
  }
  if (!(study$lor.se.var %in% colnames(data))) {
    temp <- study[["or.ucl.var", exact=TRUE]]
    if (!is.null(temp)) {
      temp  <- study[["or.conf", exact=TRUE]]
      temp  <- temp + (1-temp)/2
      zcrit <- qnorm(temp, lower.tail=TRUE)

      temp  <- log(data[, study$or.ucl.var]) - data[, study$lor.var]
      data[, study$lor.se.var] <- temp/zcrit
    }
  }

  # Compute LOR se
  if (!(study$lor.se.var %in% colnames(data))) {
    data[, study$lor.se.var] <- data[, study$or.se.var]/data[, study$or.var]
  }

  # Compute OR se
  if (!(study$or.se.var %in% colnames(data))) {
    data[, study$or.se.var] <- data[, study$lor.se.var]*data[, study$or.var]
  }

  # Compute tests and p-values
  temp <- data[, study$lor.var]/data[, study$lor.se.var]
  data[, study$test.var] <- temp
  data[, study$pvalue.var] <- 2*pnorm(abs(temp), lower.tail=FALSE)

  study$data <- data
  study

} # END: meta.setup

# Function to determine the global major and minor alleles across studies.
# For each study, the allele counts or genotype counts must be specified.
globalAlleles <- function(study.list, op=NULL) {

  # study.list         List of sublist for each study. Each sublist has the names:
  #  file              No default
  #  file.type         The default is 3
  #  delimiter         The default is "\t"
  #  header            The default is 1
  #  snp.var           The default is "SNP"
  #  major.allele      Variable for major allele
  #                    The default is "MAJOR_ALLELE"
  #  minor.allele      Variable for minor allele
  #                    The default is "MINOR_ALLELE"
  #  a.major           Major allele count (optional)
  #  a.minor           Minor allele count (optional). Must be specified
  #                     if a.major is specified.
  #  cntl00            Genotype frequency count for the major homozygous 
  #                      genotype (optional)
  #  cntl01            Genotype frequency count for the heterozygous 
  #                      genotype (optional)
  #  cntl11            Genotype frequency count for the minor homozygous 
  #                      genotype (optional)
  ###############################################################################
  # op                 List with names
  #  outfile
  #  snps              Character vector of snps names to include
  #                    or a list to scan
  #  snps.file.list    0 or 1  Set to 1 if snps is a file list

  nstudy <- length(study.list)
  op <- default.list(op, c("snps.file.list"), list(0))

  # Check each list
  for (i in 1:nstudy) {
    tlist <- default.list(study.list[[i]],
     c("file", "file.type", "header", "delimiter", "snp.var", 
       "major.allele", "minor.allele", "miss.allele"),
     list("ERROR", 3, 1, "\t", "SNP", "MAJOR_ALLELE", "MINOR_ALLELE",
          c(NA, "", " ", "-")), 
     error=c(1, 0, 0, 0, 0, 0, 0, 0))

    aflag <- (!is.null(tlist[["a.major", exact=TRUE]])) &&
             (!is.null(tlist[["a.minor", exact=TRUE]]))
    tlist$aflag <- aflag
   
    temp <- c(tlist$snp.var, tlist$major.allele, tlist$minor.allele)
    if (aflag) {
      temp <- c(temp, tlist$a.major, tlist$a.minor)
    } else {
      temp <- c(temp, tlist$cntl00, tlist$cntl01, tlist$cntl11)
    }
    tlist$keepVars <- unique(temp)

    # Check that the variables exist on the data set
    temp <- checkVars(tlist, tlist$keepVars)

    study.list[[i]] <- tlist
  } # END: for (i in 1:nstudy)

  snps <- op[["snps", exact=TRUE]]
  temp <- op[["snps.file.list", exact=TRUE]]
  if (temp) {
    snps    <- scan(snps$file, what="character", sep="\n")
    op$snps <- snps
  }

  # Read in the data, and then subset by snps and vars
  print("Reading in the data")
  for (i in 1:nstudy) {
    tlist <- study.list[[i]]
    print(tlist$name)
    temp  <- data.frame(loadData(tlist$file, tlist), stringsAsFactors=FALSE)
    temp  <- removeOrKeepCols(temp, tlist$keepVars, which=1)
    temp  <- unfactor.all(temp)
    if (!is.null(snps)) {
      ids   <- temp[, tlist$snp.var] %in% snps
      temp  <- removeOrKeepRows(temp, ids, which=1)
    }
    tlist$data <- temp

    # Check for missing snps
    if (!is.null(op[["snps", exact=TRUE]])) {
      temp <- is.na(match(op$snps, temp[, tlist$snp.var]))
      if (any(temp)) {
        print(op$snps[temp])
        temp <- paste("The above SNPs were not found for study ", i, sep="")
        print(temp)
      }
    }

    # Get the intersection of the snps for the analysis
    if (i == 1) {
      if (is.null(snps)) {
        snps <- tlist$data[, tlist$snp.var]
      } else {
        snps <- intersect(snps, tlist$data[, tlist$snp.var])
      }
    } else {
      snps <- intersect(snps, tlist$data[, tlist$snp.var])
    }

    study.list[[i]] <- tlist
  }

  # Set up and order each data frame
  for (i in 1:nstudy) {
    tlist <- study.list[[i]]
    temp  <- tlist$data
    rownames(temp) <- temp[, tlist$snp.var]
    temp  <- removeOrKeepRows(temp, snps, which=1)
    tlist$data <- temp
    study.list[[i]] <- tlist
  }

  rm(snps)
  temp <- gc()

  # Check the alleles
  print("Checking alleles")
  remove <- NULL
  for (i in 1:nstudy) {
    tlist <- study.list[[i]]
    data  <- tlist$data
    
    # Remove snps that have a missing major allele
    temp  <- data[, tlist$major.allele] %in% tlist$miss.allele
    if (any(temp)) remove <- c(remove, data[temp, tlist$snp.var])

    # If minor allele is missing, then set it to the major allele
    temp  <- data[, tlist$minor.allele] %in% tlist$miss.allele
    if (any(temp)) data[temp, tlist$minor.allele] <- data[temp, tlist$major.allele]

    # Make sure that all studies have the same unordered alleles
    a1 <- data[, tlist$major.allele] 
    a2 <- data[, tlist$minor.allele]
    if (i == 1) {
      major     <- a1
      minor     <- a2
      oneAllele <- major == minor
      n         <- length(major)
    } else {
      temp1 <- (a1 == major) | (a1 == minor)
      temp2 <- (a2 == major) | (a2 == minor)
      temp  <- !(temp1 & temp2)
      if (any(temp)) {
        ids <- (1:n)[temp]
        # Check for 1 allele
        rows <- (temp1 | temp2) & oneAllele
        rows <- (1:n)[rows]
        rows <- ids %in% rows
        if (any(rows)) {
          ids <- ids[rows]
          temp[ids] <- FALSE
        }
        rm(ids, rows)
      } 

      if (any(temp)) {
        temp1 <- paste("SNPs removed from ", tlist$name, ":", sep="")
        print(temp1)
        temp1 <- makeVector(data[temp, tlist$snp.var])
        print(temp1)
        remove <- c(remove, temp1)
      }
      rm(temp1, temp2)
      temp <- gc()
    }
    # Update
    tlist$data <- data
    study.list[[i]] <- tlist
  }

  rm(data, major, minor, a1, a2)
  temp <- gc()

  # Remove bad snps
  if (!is.null(remove)) {
    print("Removing snps")
    remove <- unique(remove)
    for (i in 1:nstudy) {
      tlist <- study.list[[i]]
      temp  <- tlist$data
      temp  <- removeOrKeepRows(temp, remove, which=-1)
      tlist$data <- temp
      study.list[[i]] <- tlist
    }
  }
  
  # Check the allele counts
  print("Checking frequency counts")
  for (i in 1:nstudy) {
    tlist <- study.list[[i]]
    data  <- tlist$data

    # Set missing values to 0
    if (tlist$aflag) {
      vars <- c(tlist$a.major, tlist$a.minor)
    } else {
      vars <- c(tlist$cntl00, tlist$cntl01, tlist$cntl11)
    }
    for (var in vars) {
      temp <- is.na(data[, var]) 
      if (any(temp)) data[temp, var] <- 0
    }
    rm(temp, var, vars)  
 
    # Update
    tlist$data <- data
    study.list[[i]] <- tlist
  }

  # Get the allele counts
  print("Getting allele counts")
  snps        <- rownames(study.list[[1]]$data)
  nsnp        <- length(snps)
  major.count <- rep.int(0, times=nsnp)
  minor.count <- rep.int(0, times=nsnp)
  for (i in 1:nstudy) {
    tlist <- study.list[[i]]
    data  <- tlist$data
    a1    <- data[, tlist$major.allele] 
    if (i == 1) {
      major <- a1
      minor <- data[, tlist$minor.allele]
    }

    # Match alleles 
    temp <- a1 == major

    if (tlist$aflag) {
      # Allele frequency counts are given
      if (any(temp)) {
        major.count[temp] <- major.count[temp] + data[temp, tlist$a.major] 
        minor.count[temp] <- minor.count[temp] + data[temp, tlist$a.minor]
      }
      if (any(!temp)) {
        major.count[!temp] <- major.count[!temp] + data[!temp, tlist$a.minor] 
        minor.count[!temp] <- minor.count[!temp] + data[!temp, tlist$a.major]
      }
    } else {
      # Genotype frequency counts are given
      if (any(temp)) {
        major.count[temp] <- major.count[temp] + 2*data[temp, tlist$cntl00] + 
                                                   data[temp, tlist$cntl01]

        minor.count[temp] <- minor.count[temp] + 2*data[temp, tlist$cntl11] + 
                                                   data[temp, tlist$cntl01]
      }
      if (any(!temp)) {
        major.count[!temp] <- major.count[!temp] + 2*data[!temp, tlist$cntl11] + 
                                                     data[!temp, tlist$cntl01]

        minor.count[!temp] <- minor.count[!temp] + 2*data[!temp, tlist$cntl00] + 
                                                     data[!temp, tlist$cntl01]
      }
    }
  }

  rm(data, a1, temp)
  temp <- gc()

  # For each snp, make sure we have both alleles
  rows <- major == minor
  if ((any(rows)) && (nstudy > 1)) {
    print("Checking minor allele")
    for (i in 2:nstudy) {
      tlist <- study.list[[i]]
      data  <- tlist$data
      a1    <- data[, tlist$minor.allele]
      temp  <- rows & (minor != a1)
      if (any(temp)) {
        minor[temp] <- a1[temp]
        rows <- major == minor
      }
      if (!any(rows)) break
      a1    <- data[, tlist$major.allele]
      temp  <- rows & (minor != a1)
      if (any(temp)) {
        minor[temp] <- a1[temp]
        rows <- major == minor
      }
      if (!any(rows)) break
    }
    rm(data, a1, rows)
  }
  
  rm(tlist, study.list)
  temp <- gc()

  # Determine the global major and minor alleles
  # Create return object
  global.major   <- rep(" ", times=nsnp)
  global.minor   <- rep(" ", times=nsnp)
  global.major.n <- rep(0, times=nsnp)
  global.minor.n <- rep(0, times=nsnp)

  temp <- major.count >= minor.count
  global.major[temp]   <- major[temp]
  global.major.n[temp] <- major.count[temp]
  global.minor[temp]   <- minor[temp]
  global.minor.n[temp] <- minor.count[temp]
  global.major[!temp]   <- minor[!temp]
  global.major.n[!temp] <- minor.count[!temp]
  global.minor[!temp]   <- major[!temp]
  global.minor.n[!temp] <- major.count[!temp]

  rm(major, minor, major.count, minor.count)
  temp <- gc()

  # Set up data frame
  ret <- data.frame(snps, global.major, global.minor, 
                    global.major.n, global.minor.n, stringsAsFactors=FALSE)
  colnames(ret) <- c("SNP", "MAJOR_ALLELE", "MINOR_ALLELE", 
                     "MAJOR_COUNT", "MINOR_COUNT")
  rownames(ret) <- ret[, "SNP"]

  # Add MAF
  ret[, "MAF"] <- ret[, "MINOR_COUNT"]/(ret[, "MAJOR_COUNT"] + ret[, "MINOR_COUNT"]) 

  out <- op[["outfile", exact=TRUE]]
  if (!is.null(out)) {
    write.table(ret, file=out, row.names=FALSE, quote=FALSE, sep="\t")
  }

  ret  

} # END: globalAlleles

# Function to flip alleles. First study is the reference study
allele.flip <- function(study.list, op=NULL) {
  
  # study.list         List of sublists with names
  #  old.major
  #  old.minor
  #  new.major
  #  new.minor
  #  snp.var
  #  flip.var
  #  outfile
  ##############################################################
  # op                 List with names
  #  snps              Character vector of snps names to include
  #                    or a list to scan
  #  snps.file.list    0 or 1  Set to 1 if snps is a file list

  nstudy <- length(study.list)
  op <- default.list(op, c("snps.file.list"), list(0))

  # Check each list
  for (i in 1:nstudy) {
    tlist <- default.list(study.list[[i]],
     c("file", "file.type", "header", "delimiter", "snp.var", 
       "old.major", "old.minor", "allele.miss", "new.major", "new.minor",
       "flip.var", "outfile"),
     list("ERROR", 3, 1, "\t", "SNP", "ERROR", "ERROR",
          c(NA, "", " ", "-"), "MAJOR_ALLELE", "MINOR_ALLELE", "FLIP", "ERROR"), 
     error=c(1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1))

     temp <- c(tlist$old.major, tlist$old.minor, tlist$new.major, tlist$new.minor)
     if (length(unique(temp)) != 4) stop("Allele names are not unique")

     if (tlist$file.type %in% c(3, 6, 8)) {
      tlist$method       <- 2
      tlist$what         <- "character"
      tlist$include.row1 <- tlist$header
    }

    study.list[[i]] <- tlist
  } # END: for (i in 1:nstudy)

  snps <- op[["snps", exact=TRUE]]
  temp <- op[["snps.file.list", exact=TRUE]]
  if (temp) {
    snps    <- scan(snps$file, what="character", sep="\n")
    op$snps <- snps
  }

  # Read in the data, and get the interscetion of the snps
  print("Reading in the data")
  for (i in 1:nstudy) {
    tlist <- study.list[[i]]
    print(tlist$name)
    data <- loadData(tlist$file, tlist)

    # Check that variables exist
    temp <- c(tlist$snp.var, tlist$old.minor, tlist$old.major)
    check.vec(temp, "vars", list(checkList=colnames(data)))

    if (!is.null(snps)) {
      ids  <- data[, tlist$snp.var] %in% snps
      data <- removeOrKeepRows(data, ids, which=1)
      rm(ids)
      gc()
    }

    # Check for missing snps
    if (!is.null(op[["snps", exact=TRUE]])) {
      temp <- is.na(match(op$snps, data[, tlist$snp.var]))
      if (any(temp)) {
        print(op$snps[temp])
        temp <- paste("The above SNPs were not found for study ", i, sep="")
        print(temp)
      }
    }

    # Get the intersection of the snps for the analysis
    if (i == 1) {
      if (is.null(snps)) {
        snps <- data[, tlist$snp.var]
      } else {
        snps <- intersect(snps, data[, tlist$snp.var])
      }
    } else {
      snps <- intersect(snps, data[, tlist$snp.var])
    }
  }

  if (is.null(snps)) {
    stop("ERROR: no common SNPs")
  }

  # Now that we have the intersection of the snps, match the alleles 
  for (i in 1:nstudy) {
    tlist  <- study.list[[i]]
    majvar <- tlist$old.major
    minvar <- tlist$old.minor
    miss   <- tlist$allele.miss
 
    print(tlist$name)
    data  <- loadData(tlist$file, tlist)
    
    # Order the data 
    rownames(data) <- data[, tlist$snp.var]
    data <- removeOrKeepRows(data, snps, which=1)

    # Let data be a data frame
    data <- data.frame(data, stringsAsFactors=FALSE)
    vars <- c(tlist$snp.var, tlist$old.major, tlist$old.minor)
    for (var in vars) data[, var]  <- unfactor(data[, var])

    # Get the reference major and minor alleles
    if (i == 1) {
      major      <- makeVector(data[, majvar])
      minor      <- makeVector(data[, minvar])
      majNotMiss <- !(major %in% miss)
      minNotMiss <- !(minor %in% miss)
      next
    } 

    # See if alleles match
    temp11 <- (major == data[, majvar])
    temp12 <- (major == data[, minvar])
    temp21 <- (minor == data[, majvar])
    temp22 <- (minor == data[, minvar])    

    temp1 <- (temp11 | temp12) & majNotMiss
    temp1[!majNotMiss] <- TRUE

    temp2 <- (temp21 | temp22) & minNotMiss
    temp2[!minNotMiss] <- TRUE

    # t1 is whether nci alleles are in new alleles
    t1 <- temp1 & temp2

    temp  <- !(data[, majvar] %in% miss)
    temp1 <- (temp11 | temp21) & temp
    temp1[!temp] <- TRUE

    temp  <- !(data[, minvar] %in% miss)
    temp2 <- (temp12 | temp22) & temp
    temp2[!temp] <- TRUE

    t2 <- temp1 & temp2

    rm(temp11, temp12, temp21, temp22, temp1, temp2, temp)
    gc()

    # Add a flip column
    data[, tlist$flip.var] <- as.numeric(!(t1 | t2))
    rm(t1, t2)
    gc()

    # Define the new alleles
    newmajvar <- tlist$new.major
    newminvar <- tlist$new.minor
    data[, newmajvar] <- data[, majvar]
    data[, newminvar] <- data[, minvar]
    temp <- data[, tlist$flip.var] == 1
    # Change alleles
    data <- changeAlleles(data, c(majvar, minvar), rows=temp,
             newVars=c(newmajvar, newminvar))
 
    ##################################################################
    # Do an error check
    print("Check for errors")
    temp11 <- (major == data[, newmajvar])
    temp12 <- (major == data[, newminvar])
    temp21 <- (minor == data[, newmajvar])
    temp22 <- (minor == data[, newminvar])

    temp1 <- (temp11 | temp12) & majNotMiss
    temp1[!majNotMiss] <- TRUE

    temp2 <- (temp21 | temp22) & minNotMiss
    temp2[!minNotMiss] <- TRUE

    # t1 is whether nci alleles are in new alleles
    t1 <- temp1 & temp2

    temp  <- !(data[, majvar] %in% miss)
    temp1 <- (temp11 | temp21) & temp
    temp1[!temp] <- TRUE

    temp  <- !(data[, minvar] %in% miss)
    temp2 <- (temp12 | temp22) & temp
    temp2[!temp] <- TRUE

    t2 <- temp1 & temp2
    rm(temp11, temp12, temp21, temp22, temp1, temp2, temp)
    gc()

    temp <- !(t1 | t2)
    rm(t1, t2)
    gc()

    if (any(temp)) {
      tsnp <- snps[temp]
      print(tsnp)
      print("ERROR with alleles for the above SNPs")
      temp <- data[, tlist$snp.var] %in% tsnp
      temp[is.na(temp)] <- FALSE
      data <- removeOrKeepRows(data, temp, which=-1)
      rm(tsnp)
      gc()
    }
    #######################################################################
    write.table(data, file=tlist$outfile, sep="\t", row.names=FALSE, quote=FALSE)
 
  } # END: for (i in 1:nstudy)

} # END: allele.flip

# Function to flip LORs. 
LOR.flip <- function(study.list, allele.list, op=NULL) {
  
  # study.list         List of sublists with names
  #  allele.var
  #  snp.var
  #  flip.var
  #  lor.var           Vector or LOR variable names
  #  new.lor           Vector of new LOR variable names
  #  outfile
  ##############################################################
  # allele.list        Type file.list
  #   snp.var
  #   allele.var
  ##############################################################
  # op                 List with names
  #  snps              Character vector of snps names to include
  #                    or a list to scan
  #  snps.file.list    0 or 1  Set to 1 if snps is a file list
  #  allele.list

  # Check for global allele list
  allele.list <- default.list(allele.list, 
           c("file", "file.type", "delimiter", "header", "allele.var", "snp.var"), 
           list("ERROR", 3, "\t", 1, "ERROR", "ERROR"), 
           error=c(1, 0, 0, 0, 1, 1)) 

  nstudy <- length(study.list)
  op <- default.list(op, c("snps.file.list", "debug"), list(0, 0))

  # Check each list
  for (i in 1:nstudy) {
    tlist <- default.list(study.list[[i]],
     c("file", "file.type", "header", "delimiter", "snp.var", 
       "allele.miss", "allele.var", "flip.var", "outfile"),
     list("ERROR", 3, 1, "\t", "SNP", 
          c(NA, "", " ", "-"), "META_MINOR_ALLELE", "FLIP.LOR", "ERROR"), 
     error=c(1, 0, 0, 0, 0, 0, 0, 0, 1))

    temp <- tlist[["new.lor", exact=TRUE]]
    if (is.null(temp)) tlist$new.lor <- tlist$lor.var
    if (length(tlist$lor.var) != length(tlist$new.lor)) {
      stop("ERROR with lor.var/new.lor")
    }

    study.list[[i]] <- tlist
  } # END: for (i in 1:nstudy)

  snps <- op[["snps", exact=TRUE]]
  temp <- op[["snps.file.list", exact=TRUE]]
  if (temp) {
    snps    <- scan(snps$file, what="character", sep="\n")
    op$snps <- snps
  }
  debug <- op$debug

  # Read in the global alleles
  print("Reading global allele file")
  t00      <- proc.time()
  data     <- getColumns(allele.list, c(allele.list$snp.var, allele.list$allele.var))
  g.snp    <- data[[allele.list$snp.var]]
  g.allele <- data[[allele.list$allele.var]]
  names(g.allele) <- g.snp
    
  if (is.null(snps)) {
    snps <- g.snp
  } else {
    snps <- intersect(snps, g.snp)
  }
  rm(g.snp, data)
  temp <- gc()

  # Read in the data, and get the intersection of the snps
  if (debug) t00 <- debug.time(t00, str="Reading in the data")
  for (i in 1:nstudy) {
    tlist <- study.list[[i]]
    print(tlist$name)
    if (tlist$file.type %in% c(3, 6, 8)) {
      tlist$method       <- 2
      tlist$what         <- "character"
      tlist$include.row1 <- tlist$header
    }
    data <- loadData(tlist$file, tlist)

    # Check that variables exist
    vars <- c(tlist$snp.var, tlist$allele.var, tlist$lor.var)
    check.vec(vars, "vars", list(checkList=colnames(data)))
    rownames(data) <- data[, tlist$snp.var] 

    # Subset the data
    if (debug) t00 <- debug.time(t00, str="Subsetting the data")
    temp <- data[, tlist$snp.var] %in% snps
    data <- removeOrKeepRows(data, temp, which=1)

    # Let data be a data frame
    if (debug) t00 <- debug.time(t00, str="Create data frame")
    data <- data.frame(data, stringsAsFactors=FALSE)
    for (var in vars) data[, var] <- unfactor(data[, var]) 
    for (var in tlist$lor.var) data[, var] <- as.numeric(data[, var])

    if (debug) t00 <- debug.time(t00, str="Check for bad SNPs")
    temp <- !(data[, tlist$snp.var] %in% snps)
    if (any(temp)) {
      #data[, "NotInAlleleFile"] <- 0
      #data[temp, "NotInAlleleFile"] <- 1
      print("Some SNPs were not found in allele file. These SNPs will be removed.")
      data <- removeOrKeepRows(data, temp, which=-1)
    }

    # Initialize
    if (debug) t00 <- debug.time(t00, str="Initialize new variables")

    data[, tlist$flip.var] <- 0
    vars <- tlist$lor.var
    new  <- tlist$new.lor
    for (j in 1:length(vars)) {
      data[, new[j]] <- data[, vars[j]]
    }
 
    # Order the allele vector
    if (debug) t00 <- debug.time(t00, str="Order the allele vector")

    #temp   <- makeVector(data[, tlist$snp.var])
    #allele <- g.allele[temp]
    temp   <- match(makeVector(data[, tlist$snp.var]), names(g.allele)) 
    temp   <- temp[!is.na(temp)]
    allele <- g.allele[temp]
    temp   <- allele != data[, tlist$allele.var]
  
    if (any(temp)) {
      if (debug) t00 <- debug.time(t00, str="Flip LOR")
      data[temp, tlist$flip.var] <- 1
      vars <- tlist$lor.var
      new  <- tlist$new.lor
      for (j in 1:length(vars)) {
        data[temp, new[j]] <- -data[temp, vars[j]]
      }
    }

    if (debug) t00 <- debug.time(t00, str="Writing table to file")
    write.table(data, file=tlist$outfile, sep="\t", row.names=FALSE, quote=FALSE)
  }

  0

} # END: LOR.flip

# MiMa: An S-Plus/R function to fit meta-analytic mixed-, random-, and fixed-effects models
mima <- function(yi, vi, mods, method="REML", threshold=0.00001, maxiter=100, alpha=0.05, digits=4, fe="no", verbose="no", out="no") {

	k			<- length(yi)
	intrcpt	<- rep(1,k) 
	X			<- cbind(intrcpt, mods)
	y			<- as.matrix(yi)
	p			<- dim(X)[2] - 1				### number of moderators in the model

	### function to obtain the trace of a matrix 
	tr <- function(X) {
		sum(diag(X))
	}

	if (fe == "no") {

		if (method == "HE") {
			M		<- X %*% solve(t(X) %*% X) %*% t(X)
			RSS	<- t(y) %*% ( diag(k) - M ) %*% y
			vart	<- ( RSS - tr( (diag(k)-M) %*% diag(vi) ) ) / ( k-p-1 )
		}

		if (method == "DL") {
			wi		<- 1/vi
			W		<- diag(wi)
			b		<- solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% y
			RSS	<- t(y - X %*% b) %*% W %*% (y - X %*% b)
			vart	<- ( RSS - (k-p-1) ) / ( sum(wi) - tr( t(X) %*% W %*% W %*% X %*% solve( t(X) %*% W %*% X ) ) )
		}

		if (method == "SH") {
			vart0	<- sum( (yi - mean(yi))^2 ) / k
			wi		<- 1/(vi + vart0) 
			W		<- diag(wi)
			P		<- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
			vart	<- (vart0/(k-p-1)) * t(y)%*%P%*%y
		}

		conv		<- 1
		itprec	<- threshold
		change	<- 1
		iter		<- 0

		if (method == "REML" || method == "ML" || method == "EB") {

			vart		<- var(yi) - 1/k*sum(vi)	### initial estimate = HE estimator in RE model
			vart[vart < 0] <- 0						### set to zero in case the initial estimate is negative

			while (change > itprec) {
				if(verbose == "yes") cat("Iteration:", iter, " Estimate of (Residual) Heterogeneity:", round(vart, 8), "\n")
				iter	<- iter + 1
				varm	<- vart
				wi		<- 1/(vi + vart)
				W		<- diag(wi)
				V		<- diag(vi)
				P		<- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
				if (method == "REML") {
					adj	<- solve( tr(P%*%P) ) %*% ( t(y)%*%P%*%P%*%y - tr(P) )
				}
				if (method == "ML") {
					adj	<- solve( tr(W%*%W) ) %*% ( t(y)%*%P%*%P%*%y - tr(W) )
				}
				if (method == "EB") {
					adj	<- solve( tr(W) ) %*% ( (k/(k-p-1)) %*% t(y)%*%P%*%y - k )
				}
				while (vart + adj < 0) { 
					adj <- adj / 2
				}
				vart		<- vart + adj
				change	<- abs(varm - vart)
				if (iter > maxiter) {
					conv    <- 0
					break
				}
			}
		}

		if (conv == 0) {
                  return(NULL)
			cat("Fisher scoring algorithm did not converge\nTry increasing maxiter or use a different estimation method\n")
			break
		}

	}

	if (fe == "yes") vart <- 0

	wi		<- 1/vi
	W		<- diag(wi)
	b		<- solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% y
	vb		<- solve(t(X) %*% W %*% X)
	QE		<- sum(wi * yi^2) - t(b) %*% solve(vb) %*% b
	QEp	<- 1-pchisq(QE, df=k-p-1)

	vart[vart < 0] <- 0
	wi		<- 1/(vi + vart)
	W		<- diag(wi)
	b		<- solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% y
	vb		<- solve(t(X) %*% W %*% X)

	if (p > 0) {
		QME	<- t(b)[2:(p+1)] %*% solve(vb[2:(p+1),2:(p+1)]) %*% b[2:(p+1)]
		QMEp	<- 1 - pchisq(QME, df=p)
	} else {
		QME	<- NA
		QMEp	<- NA
	}

	zvals	<- b/sqrt(diag(vb))
	zp		<- 2*(1-pnorm(abs(zvals)))
	lbsci	<- b - qnorm(1-alpha/2) * sqrt(diag(vb))
	ubsci	<- b + qnorm(1-alpha/2) * sqrt(diag(vb))

	if (out == "yes") {
            heterogeneity.test <- list(QE, k-p-1, QEp) 
            names(heterogeneity.test)<-c("Het_test", "df", "p_het")
            effect.estimate <- data.frame(cbind(b, sqrt(diag(vb)), zvals, zp))
		names(effect.estimate) <- c("estimate", "SE", "zval", "pval")
            outlist<-list(heterogeneity.test=heterogeneity.test,
				effect.estimate=effect.estimate)

	##	outlist <- list(vart, b, vb)
	##	names(outlist) <- c("vart", "b", "vb")
		outlist
	} else {

	cat("\n")
	cat("Estimate of (Residual) Heterogeneity:", round(vart,digits))
	cat("\n\n")

	cat("Test for (Residual) Heterogeneity:")
	cat("\n\n")
	cat("QE      = ", round(QE,digits))
	cat("\n")
	cat("df      = ", k-p-1)
	cat("\n")
	cat("p-value = ", round(QEp,digits))
	cat("\n\n")

	cat("Parameter Estimates:")
	cat("\n\n")
	print(round(b,digits))
	cat("\n")

	cat("Variance-Covariance Matrix of Parameter Estimates:")
	cat("\n\n")
	print(round(vb,digits))
	cat("\n")

	cat("Omnibus Test of all Moderators:")
	cat("\n\n")
	cat("QME     = ", round(QME,digits))
	cat("\n")
	cat("df      = ", p)
	cat("\n")
	cat("p-value = ", round(QMEp,digits))
	cat("\n\n")

	cat("Individual Moderator Tests:")
	cat("\n\n")
	ztest <- cbind(b, sqrt(diag(vb)), zvals, zp, lbsci, ubsci)
	ztest	<- data.frame(ztest)
	names(ztest) <- c("estimate", "SE", "zval", "pval", "CI_L", "CI_U")
	print(round(ztest,digits))
	cat("\n")

	}
} # END mima

# Function for type 3 meta-analysis (mima.R)
meta.mima.scan <- function(file.list, op=NULL) {

  # file.list
  #   lor.var   Character vector of variable names for study specific LORs
  #             No default
  #   se.var    Character vector of variable names for study specific SEs
  #             No default
  #   snp.var   No default
  #########################################################################
  # op          List
  #  outfile
  #  rows       Vector of row numbers to use, or character vector of
  #             SNP names.
  #             The default is NULL.
  #  fe         "no" or "yes" for fixed effects analysis
  #             The default is "no"
  #  mods       Vector of moderators
  #             The default is c()
  #########################################################################

  file.list <- default.list(file.list, 
               c("file", "file.type", "delimiter", "header", 
                 "lor.var", "se.var", "snp.var"), 
               list("ERROR", 3, "\t", 1, "ERROR", "ERROR", "ERROR"), 
               error=c(1, 0, 0, 0, 1, 1, 1))

  op <- default.list(op, c("fe", "mods"), list("no", c()), error=c(0, 0))
  op$fe <- tolower(op$fe)
  if (op$fe != "yes") op$fe <- "no"

  # Read in the data
  print("Reading the data")
  x <- read.table(file.list$file, sep=file.list$delimiter, header=file.list$header)
  
  print("Unfactoring columns")
  x <- unfactor.all(x)
  lorVars <- file.list$lor.var
  seVars  <- file.list$se.var
  snpVar  <- file.list$snp.var

  # Subset the data
  rows <- op[["rows", exact=TRUE]]
  if (!is.null(rows)) {
    print("Subsetting data")
    if (is.numeric(rows)) {
      x <- removeOrKeepRows(x, rows, which=1)
    } else {
      # Assuming SNP names
      temp <- x[, snpVar] %in% rows  
      x <- removeOrKeepRows(x, temp, which=1)
    }
  }
  rm(rows)
  gc()
  
  for (var in c(lorVars, seVars)) x[, var] <- as.numeric(x[, var])

  # Get the snp names
  snps <- makeVector(x[, snpVar])

  # Set up return data frame
  nr   <- nrow(x)
  cols <- c("C", "N", "N", "N", "N", "N", "N")
  cnames <- c("SNP", "Heter.test", "Heter.df", "Heter.pvalue", 
              "Estimate", "SE", "Pvalue") 
  ret <- initDataFrame(nr, cols, rownames=NULL, colnames=cnames,
                  initChar="", initNum=NA)
  ret[, 1] <- snps

  x <- removeOrKeepCols(x, snpVar, which=-1)

  rm(snps, snpVar)
  temp <- gc()
  
  # Moderators
  mods <- op$mods
  fe   <- op$fe

  # Loop over each row
  print("Begin analysis")
  for (i in 1:nr) {
    lor <- makeVector(x[i, lorVars])
    se  <- makeVector(x[i, seVars])

    # Get the variances
    se <- se*se

    # Call mima
    out <- try(mima(lor, se, mods, out="yes", fe=fe), silent=TRUE)

    if (class(out) == "try-error") next
    if (is.null(out)) next

    temp <- out$heterogeneity.test 
    ret[i, "Heter.test"]   <- temp$Het_test
    ret[i, "Heter.df"]     <- temp$df
    ret[i, "Heter.pvalue"] <- temp$p_het
    temp <- out$effect.estimate 
    ret[i, "Estimate"]     <- temp[1, "estimate"]
    ret[i, "SE"]           <- temp[1, "SE"]
    ret[i, "Pvalue"]       <- temp[1, "pval"]
  }

  # Output
  temp <- op[["outfile", exact=TRUE]]
  if (!is.null(temp)) writeTable(ret, temp)

  ret

} # END: meta.mima.scan

# Function for fixed effects meta analysis, combining allele.flip, lor.flip, and meta.analysis
meta <- function(study.list, op=NULL) {

  # study.list        List of sublists with names
  #  file                 No default
  #  major.var            Name of major allele. No default
  #  minor.var            Name of minor allele. No default
  #  snp.var              Name of snp variable. No default
  #  lor.var              Name of log-odds ratio variable. No default
  #  lor.se.var           Name of the variable for the standard error of the log-odds ratio. No default
  # op                List with names
  #  outfile              Default is NULL
  #  temp.list
  #  method               Default is "1a"
  #  allele.list
  #  copy.vars
  #  copy2.vars
  #  n.case               (For method = 2)
  #  n.control            (For method = 2)
  #  add.vars             Variables to appear in outfile for each study.
  #                       ex. c("lor.var", "lor.se.var", "or.var")
  #                       The default is NULL
  #  add.names            Variable names for add.vars.
  #                       ex. c("LOR", "LOR.SE", "OR")
  #                       The default is NULL

  temp.list <- op[["temp.list", exact=TRUE]]
  temp.list <- check.temp.list(temp.list)
  dir       <- checkForSep(temp.list$dir)
  id        <- temp.list$id
  delete    <- temp.list$delete
  rm(temp.list)
  gc()

  op <- default.list(op, c("method"), list("1a"))
  nstudy <- length(study.list)
  all.files <- NULL

  # Check each list
  print("########################")
  print("Checking the input lists")
  print("########################")
  for (i in 1:nstudy) {
    tlist <- study.list[[i]]
    tlist <- default.list(tlist,
     c("file", "snp.var", "major.var", "minor.var", "name", "lor.var"),
     list("ERROR", "ERROR", "ERROR", "ERROR", paste("study_", i, sep=""), "ERROR"), 
     error=c(1, 1, 1, 1, 0, 1))

    vars <- c(tlist$snp.var, tlist$major.var, tlist$minor.var)
    tlist <- check.file.list(tlist, op=list(vars=vars))
  
    # Add options
    tlist$old.major <- tlist$major.var
    tlist$old.minor <- tlist$minor.var
    tlist$new.major <- paste(tlist$major, "_new", sep="")
    tlist$new.minor <- paste(tlist$minor, "_new", sep="")
    tlist$flip.lor  <- 0

    # temporary file
    temp <- paste(tlist$name, "_", id, "_a", sep="")
    tmp <- getTempfile(dir, prefix=temp, ext=".txt.xls")
    tlist$outfile <- tmp
    all.files <- c(all.files, tmp)

    if (i == 1) {
      tlist$allele.var <- tlist$minor.var
    } else {
      tlist$allele.var <- tlist$new.minor
    }
    tlist$new.lor <- paste(tlist$lor.var, "_new", sep="")

    study.list[[i]] <- tlist

  } # END: for (i in 1:nstudy)

  # Check for allele list. If not specified, use the first study
  allele.list <- op[["allele.list", exact=TRUE]]
  aflag <- !is.null(allele.list)
  if (aflag) {
    allele.list <- default.list(allele.list, c("snp.var", "allele.var"), 
                    list("ERROR", "ERROR"), error=c(1, 1))
    vars <- c(allele.list$snp.var, allele.list$allele.var)
    allele.list <- check.file.list(allele.list, op=list(vars=vars))
  } else {
    # Use the first study
    allele.list <- study.list[[1]]
  }
  
  # Flip alleles
  print("########################")
  print("Calling allele.flip")
  print("########################")
  ret <- allele.flip(study.list, op=op) 
 
  # Define the new files
  for (i in 2:nstudy) {
    tlist <- study.list[[i]]
    tlist$file      <- tlist$outfile
    tlist$file.type <- 3
    tlist$delimiter <- "\t"
    tlist$header    <- 1
    temp <- paste(tlist$name, "_", id, "_lor", sep="")
    tmp <- getTempfile(dir, prefix=temp, ext=".txt.xls")
    tlist$outfile <- tmp
    all.files <- c(all.files, tmp)
    study.list[[i]] <- tlist
  }
  
  # Flip lor
  print("########################")
  print("Calling LOR.flip")
  print("########################")
  ret <- LOR.flip(study.list, allele.list, op=op)

  op$allele.list <- NULL
  rm(allele.list, temp, tmp)
  gc()

  # Define the new files
  for (i in 1:nstudy) {
    tlist <- study.list[[i]]
    tlist$file      <- tlist$outfile
    tlist$file.type <- 3
    tlist$delimiter <- "\t"
    tlist$header    <- 1
    tlist$lor.var   <- tlist$new.lor
    tlist$new.lor   <- NULL
    study.list[[i]] <- tlist
  }

  # Run meta-analysis
  print("########################")
  print("Calling meta.analysis")
  print("########################")
  ret <- meta.analysis(study.list, op=op)

  if (delete) {
    for (tmp in all.files) file.remove(tmp)
  }

  ret

} # END: meta

# Function to read in  obs file and convert to a data frame
convert.obsfile <- function(obs.outfile, op=NULL) {

  op <- default.list(op, c("add.CI", "alpha", "gene.list", "ci.digits"), 
                     list(1, 0.05, NULL, 4))

  fid    <- file(obs.outfile, "r")
  x1     <- scan(fid, what="character", sep=",", nlines=1)
  x2     <- scan(fid, what=double(0), sep=",", nlines=1)
  x3     <- scan(fid, what=double(0), sep=",", nlines=1)
  x4     <- scan(fid, what=double(0), sep=",", nlines=1)
  x5     <- scan(fid, what=double(0), sep=",", nlines=1)
  y.cont <- scan(fid, what=double(0), sep=",", nlines=1)
  close(fid)
  x <- cbind(x1, x2, x3, x4, x5)
  rm(x1, x2, x3, x4, x5)
  gc()

  x <- as.data.frame(x, stringsAsFactors=FALSE)
  colnames(x) <- c("SNP", "Pvalue", "Pvalue.flag", "Beta", "SE")

  if (op$add.CI) {
    zcrit <- qnorm(1 - (op$alpha)/2)
    beta  <- as.numeric(x[, "Beta"])
    se    <- as.numeric(x[, "SE"])
    if (!y.cont) x[, "OR"] <- exp(beta)
    l <- beta - 1.96*se
    u <- beta + 1.96*se
    if (!y.cont) {
      l <- exp(l)
      u <- exp(u)
    }
    l <- round(l, digits=op$ci.digits)
    u <- round(u, digits=op$ci.digits)
    x[, "CI"] <- paste("(", l, ", ", u, ")", sep="")
  } 

  # Append genes
  gene.list <- op[["gene.list", exact=TRUE]]
  if (!is.null(gene.list)) {
    gene.list <- default.list(gene.list, c("snp.var", "gene.var", "chrm.var"), 
                    list("SNP", "Gene", "Chr"))
    vars <- c(gene.list[["gene.var", exact=TRUE]], gene.list[["snp.var", exact=TRUE]])
    gene.list <- check.file.list(gene.list, op=list(exist=1, vars=vars))
    gene.list$return.cols <- 1
    cols <- getNcols(gene.list$file, gene.list)
    if (gene.list$chrm.var %in% cols) {
      chr <- gene.list$chrm.var
    } else {
      chr <- NULL
    }
    gene.list$vars   <- c(gene.list$gene.var, chr)
    gene.list$id.var <- gene.list$snp.var
    x <- addColumn(x, "SNP", gene.list)
  }

  x

} # END: convert.obsfile
