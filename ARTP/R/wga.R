# History Feb 14 2008  Initial coding
#         Feb 19 2008  Add function convert.GenABEL
#                      Change recode.vec to recode.geno
#         Feb 22 2008  Add the name "file" to snp.list and remove
#                      the names "prefix" and "extension". 
#                      Change "in.dir" to "dir" in snp.list
#         Mar 04 2008  Make convert.snpMatrix more efficient
#         Mar 07 2008  Allow snps from different chromosomes to be part
#                      of the same snpMatrix data object.
#                      Move chrms, start.vec, stop.vec into snp.list
#         Mar 12 2008  Add seperate functions to check each list
#                      Add snpNames field to snp.list
#         Mar 13 2008  Change the purpose of the chrms vector in snp.list
#                      Allow reading other flat files
#         Mar 14 2008  Change dir and file fields in locusMap.list
#                      Add snp.list$row1.omitN
#                      Add funtion readPhenoData
#                      Change the function getPhenoData
#                      Add readTable function.
#         Mar 17 2008  Change recode.geno
#                      Rename readTable to table2Format
#                      Rename readPhenoData to readTable
#                      Allow phenotype data and locus map data to be type 1
#         Mar 18 2008  Check for numeric vectors when there is no header.
#                      Do not always check the phenotype list for snpMatrix
#                      Add function to transform a sas file (SAS2Format).
#         Mar 19 2008  Add getLocusMap to read the locus map data set
#                      In convert.GenABEL, delete temporary files as they
#                      are no longer needed instead of at the end.
#         Mar 20 2008  Make scanFile more general.
#                      Add file.type = 4 (SAS data sets)
#         Mar 21 2008  Change in convert.GenABEL when adding the snp data.
#                      Read in the temporary files as lines instead of each
#                      element at a time.
#                      Allow readTable not to have an id variable
#         Mar 24 2008  Use loadData function to load the phenotype and
#                      locusMap data.
#                      Allow an option to pass in the snp data file into
#                      convert.snpMatrix
#                      Change in convert.GenABEL to subset by subject ids.
#                      Change in convert.snpMatrix to subset by subject ids.
#         Mar 25 2008  Add createDummy function and createDummy option in the
#                      phenotype list.
#         Mar 26 2008  Add option to save the snp output files
#                      Re-write sas2Format
#         Mar 27 2008  Get the order of the snps in the snp files
#                      Break functions up ino seperate files
#         Mar 31 2008  Change the snpNames option
#         Apr 01 2008  Add error checks.
#         Apr 03 2008  Change in getPhenoData to check if the subject ids
#                      are unique. If not, then create a variable with
#                      unique ids. Return a list of the data and new var.
#                      Change is complete for GenABEL.
#         Apr 04 2008  Make changes in convert.snpMatrix for non-unique
#                      subject ids.
#                      Change is complete for snpMatrix
#         Apr 05 2008  Change in getDataObject for which = 1.
#                      Change the subject ids when the ids are not unique
#         Apr 08 2008  Change in loadData.type1 to return the subject ids
#                      in the order they appear in the data
#         Apr 09 2008  Add getSubjIds function
#                      Set pdata.ids to be the subject ids even when the
#                      ids are unique in getPheno.info
#         Apr 11 2008  Add getDelimiter
#         Apr 15 2008  Change to sas.list
#         Apr 17 2008  Put the names delete, id, and dir in sas.list
#                      Change in getPhenoData
#                      Put function addToPhenoData in wga_unused.R
#         Apr 19 2008  Change update.snpNames function.
#         Apr 30 2008  Add option to remove rows with missing values
#                      in getPhenoData.
#         May 02 2008  Check for same number of genotypes for each SNP.
#         May 08 2008  Add option so that loadData.type1 returns the
#                      phenotype data.
#                      Remove snp.list$snpNames after the snpNames
#                      vector is defined.
#         May 09 2008  Add inheritance option
#         May 12 2008  Remove remove.miss from getPhenoData
#         Jun 10 2008  Change inheritance to genetic.model
#         Jun 27 2008  Have getLocusMap call getColumns
#         Jun 27 2008  Call update.snp.list to get all snp ids to use
#         Jun 28 2008  Extract snpNames after loadData is called
#                      Add loadData.stream
#         Jul 08 2008  When subjects are not found, save them in a 
#                      global object bad.subject.ids.global
#         Jul 30 2008  Add function to get the intersection of the
#                      subject ids on the phenotype and genotype files.
#         Aug 14 2008  Move intersectSubIds call into getData.1
#         Sep 17 2008  Fix bug in intersectSubIds for streamed input
#         Oct 02 2008  Set snpNames and snpNames.list to NULL in 
#                      intersectSubIds function
#         Oct 10 2008  Add code not to recode genotypes
#         Oct 29 2008  Recode before subsetting the subjects (changes in
#                      loadData.type1 and orderSNP
#         Dec 11 2008  Add code for formulas in pheno.list
#         Dec 15 2008  Move snpMatrix, GenABEL code to wga_unused.R
#         Feb 06 2009  Change in getPheno.info to return the subject ids
#                      which are controls
#         Feb 19 2009  Add MAF calculation in loadData.type1
#         Mar 05 2009  Do not produce an error if cc.var is not specified in
#                      loadData.type1
#                      Compute MAF in orderSNP
#         Mar 20 2009  Check cc.var is 0-1 in getPheno.info
#         Apr 07 2009  Leave ids as variable in the phenotype data
#         May 01 2009  Update getLocusMap function to return snps from a 
#                      given chromosome in a selected range.
#         Jul 21 2009  Call checkDelimiter function in getData.1
#         Aug 10 2009  Do not stop if snpNames were not found. Add error option.
#                      Use options out.miss and out.delimiter in snp.list
#         Sep 17 2009  Add option in getDelimiter for the output snp data
#         Oct 16 2009  Set up changes:
#                      getPhenoData
#                      loadData.type1
#                      getData.1
#                      getPheno.info
#                      intersectSubIds
#         Oct 16 2009  Change in recode.geno
#         Dec 30 2009  Add code for getting variables from formulas
#         Feb 25 2010  Add genotype frequencies and number of missing genotypes
#                      to loadData.type1
#         Mar 25 2010  Add gene.var to getLocusMap
#         Apr 08 2010  Update loadData.type1 for file.type = 9, 10
#         Apr 12 2010  Change getSNPLoadList for file.type = 9, 10
#         Apr 27 2010  Add code for duplicated pheno ids
#         Jun 09 2010  Return imputed vector for imputed data in loadData.type1
#         Oct 18 2010  With an empty space missing value ("") check the last genotype. If missing,
#                      then a missing value needs to be added.

# Function to read the locus map data set
# This function return a list with the names "snp", "chrm", and "loc".
# "snp" and "chrm" are character vectors, and "loc" is a numeric vector.
getLocusMap <- function(file, locusMap.list, temp.list=NULL, op=NULL) {

  # file           File to read
  #                No default
  # locusMap.list  List of options (See locusMap.list.wordpad)
  #                No default
  # op             List with names chrm, start, stop

  # Set up the input arguments to getColumns()
  vars <- c(locusMap.list$snp.var, locusMap.list$chrm.var,
            locusMap.list$loc.var, locusMap.list$alleles.var,
            locusMap.list$gene.var)
  file.list      <- locusMap.list
  file.list$file <- file

  # Check for errors
  chrmFlag <- !is.null(locusMap.list[["chrm.var", exact=TRUE]])
  locFlag  <- !is.null(locusMap.list[["loc.var", exact=TRUE]])
  snpFlag  <- !is.null(locusMap.list[["snp.var", exact=TRUE]])
  allFlag  <- !is.null(locusMap.list[["alleles.var", exact=TRUE]])
  geneFlag <- !is.null(locusMap.list[["gene.var", exact=TRUE]])

  opFlag   <- !is.null(op)

  data <- getColumns(file.list, vars, temp.list=temp.list)

  # Return a list
  ret <- list() 
  
  if (snpFlag) ret$snp     <- data[[locusMap.list$snp.var]]
  if (chrmFlag) ret$chrm   <- data[[locusMap.list$chrm.var]]
  if (locFlag) ret$loc     <- as.numeric(data[[locusMap.list$loc.var]])
  if (allFlag) ret$alleles <- data[[locusMap.list$alleles.var]]
  if (geneFlag) ret$gene   <- data[[locusMap.list$gene.var]]


  # Determine if only certain snps in a range is desired
  if (opFlag) {
    temp <- rep(TRUE, times=length(data[[1]]))
    chr <- getListName(op, "chrm")
    if ((!is.null(chr)) & (chrmFlag)) {
      temp <- temp & (ret$chrm %in% chr)
    }
    start <- getListName(op, "start")
    if ((!is.null(start)) & (locFlag)) {
      temp <- temp & (ret$loc >= as.numeric(start))
    }
    stop <- getListName(op, "stop")
    if ((!is.null(stop)) & (locFlag)) {
      temp <- temp & (ret$loc <= as.numeric(stop))
    }
    temp[is.na(temp)] <- FALSE
    if (snpFlag) ret$snp     <- ret$snp[temp]
    if (chrmFlag) ret$chrm   <- ret$chrm[temp]
    if (locFlag) ret$loc     <- ret$loc[temp]
    if (allFlag) ret$alleles <- ret$alleles[temp]
    if (geneFlag) ret$gene   <- ret$gene[temp]
  }

  ret
 
} # END: getLocusMap

# Function to read and modify the phenotype data based on which
#  package is being used. 
getPhenoData <- function(p, temp.list=NULL) {

  # The options sex.var, sex.female, and sex.male are used
  #  with which = 2 and 3 (see below).
  
  # p is a list with the folowing fields:
  #############################################################
  # file           Phenotype data file. This file must have an
  #                id variable.
  #                No default
  # file.type      1, 3, 4  1 is for an R object file created with the
  #                save() function. 3 is for a table that will be read in
  #                with read.table(). 4 is for a SAS data set.
  #                The default is 3
  # header         0 or 1 . Set to 0 if the file does not contain a header.
  #                The default is 1.
  # delimiter      The delimiter in the table.
  #                The default is ""
  # id.var         Name or column number of the id variable. 
  #                No default.
  # keep.vars      Vector of variable names or column numbers to keep. 
  #                The default is NULL, so that all variables will be kept.
  # remove.vars    Vector of variable names or column numbers to remove. 
  #                The default is NULL, so that all variables will be kept.
  #                Both use.vars and remove.vars cannot be specified. 
  # factor.vars    Vector of variable names or column numbers to convert
  #                into factors.
  #                The default is NULL.
  # make.dummy     0 or 1 to make dummy variables for factor.vars
  #                The default is 0.
  # keep.ids       Vector of ids or row numbers to keep.
  #                The default is NULL.
  # remove.ids     Vector of ids or row numbers to remove.
  #                Both keep.ids and remove.ids cannot be specified
  #                The default is NULL.
  # in.miss        Vector of character strings to define the missing values.
  # remove.miss    0 or 1 to remove rows with missing values
  #                The default is 0.
  # sas.list       See below
  ###############################################################

  # Check the list
  p <- check.pheno.list(p)

  id.var <- p$id.var

  if (!is.null(p$keep.vars)) {
    if (!(p$id.var %in% p$keep.vars)) stop("ERROR: keep.vars must contain id.var")
  } 

  # Set options 
  p$transpose    <- 0
  p$start.row    <- 1
  p$stop.row     <- -1
  p$stream       <- 0
  p$include.row1 <- p$header
  p$snpNames     <- NULL
  p$method       <- 1
  if (p$file.type == 4) {
    p$sas.list$delimiter <- "|"
    temp.list <- check.temp.list(temp.list)
    p$sas.list$temp.list <- temp.list
  }

  # Check if missing values are to be removed.
  # NOTE: a formula can create missing values when the variable
  #  does not have any: log(x)
  removeFlag <- p$remove.miss

  # Check if any formulas are to be applied
  formulas <- getFormulas(p)
  formFlag <- length(formulas) && removeFlag

  if (formFlag) p$remove.miss <- 0

  # The data has been loaded, call readTable
  x <- readTable(p$file, p) 
  p$data <- NULL

  if (formFlag) {
    # Update the data for any formulas
    x <- applyFormulas(x, formulas)

    # Remove missing values
    vars <- getAllVars(p)
    x    <- removeMiss.vars(x, vars=vars, miss=p$in.miss)
    
  } else {
    if (removeFlag) x <- removeMiss(x, miss=p$in.miss)
  }

  nr <- nrow(x)

  # Get the id var
  id <- x[, p$id.var]

  # Check if the ids are unique
  orig.id <- NULL
  uniq.id <- unique(id)
  n.uid   <- length(uniq.id)
  if (n.uid != nr) {
    # Print a warning message
    warning("The subject ids are not all unique")
  }
  if (p$unique.ids) {
    
    orig.id <- p$orig.id.var
    
    x[, orig.id] <- x[, id.var]
 
    # Get the new ids
    count        <- rep(0, n.uid)
    names(count) <- uniq.id
    temp         <- x[, orig.id]
    dups         <- duplicated(x[, id.var])
    if (is.factor(temp)) temp <- as.character(levels(temp))[temp]
    for (i in 2:nr) {
      if (dups[i]) {
        name    <- temp[i]
        cntname <- count[name]
        if (cntname) {
          temp[i]     <- paste(name, "_", cntname, sep="")
          count[name] <- cntname + 1
        } else {
          count[name] <- 2
        }
      }
    }

    x[, id.var] <- temp
  } # END: if (p$unique.ids)

  # Set the row names
  #rownames(x) <- as.character(x[, id.var])
 
  # Return a list
  list(data=x, orig.id=orig.id)

} # END: getPhenoData

# Function to return the snp data for which = 1
loadData.type1 <- function(snp.list, pheno.list, temp.list, op=NULL) {

  # snp.list
  # pheno.list
  #####################################################
  # op               A list with the following names
  #  include.row1    0 or 1 to include the header
  #                  The default is 1
  #  include.snps    0 or 1 to include the snp names as the
  #                  first element of field.
  #                  The default is 0.
  #  return.type     1 or 2   1 is a vector of characters
  #                  2 is a matrix. The returned matrix will have
  #                  the column names as subject ids and row names
  #                  as snpnames.
  #                  The default is 1
  #  missing         0 or 1 to return a logical vector to determine
  #                  which snps had missing values
  #                  The default is 1.
  #  snpNames        0 or 1 to return a vector of snp names
  #                  The default is 1.
  #  subjIds         0 or 1 to return a vector of subject ids corresponding
  #                  to the order they appear. (Subjects are columns)
  #                  The default is 0
  #  orderByPheno    0 or 1 to order the columns of the snps according
  #                  to the order of the subject ids in the phenotype
  #                  data. 
  #                  The default is 1.
  #  return.pheno    0 or 1 to return the phenotype data
  #                  The default is 1.
  #  useControls     0 or 1 to use only the subset of controls to determine
  #                  the minor allele. If set to 1, then there must be a name
  #                  "cc.var" in pheno.list which gives the name of the 0-1
  #                  case-control variable.
  #                  The default is 1.
  #  MAF             0 or 1 to compute the MAF for each SNP. If useControls = 1,
  #                  then only the controls will be used.
  #                  The default is 0.
  #  alleles         0 or 1 to return the major/minor alleles
  #                  The default is 1
  #  genoFreqs       0 or 1 to return genotype frequencies as a string 0-1-2
  #                  The default is 0
  #  n.miss          0 or 1 to return number of missing genotypes
  #                  The default is 0.

  # Check the lists. pheno.list is checked in getPheno.info
  snp.list  <- check.snp.list(snp.list)
  temp.list <- check.temp.list(temp.list)

  # Check the options list
  if (is.null(op)) op <- list()
  op <- default.list(op, 
         c("include.row1", "include.snps", "return.type", "missing",
            "snpNames", "subjIds", "orderByPheno", "return.pheno", 
            "useControls", "MAF", "stopOnError", "alleles",
            "genoFreqs", "n.miss"),
         list(1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0))
  if (op$return.type == 2) op$include.row1 <- 1

  # Check for error condition
  if (op$useControls) {
    ccvar <- getListName(pheno.list, "cc.var")
    if (is.null(ccvar)) {
      #print("##########################################################")
      #print("NOTE: pheno.list$cc.var not specified")
      #print("##########################################################")
      op$useControls <- 0
    }
  }

  # Character vector containing the snp data. First element contains the
  # subject ids. Remaining elements contain the snp name and genotypes.
  # This character vector must be delimited.
  i              <- 1
  data.obj       <- NULL
  delimiter      <- getDelimiter(snp.list)
  in.miss        <- get.in.miss(snp.list)
  heter.codes    <- snp.list$heter.codes 
  mode           <- snp.list$genetic.model
  codes          <- getInheritanceVec(mode, recode=snp.list$recode)
  missFlag       <- op$missing
  phenoData.list <- NULL

  # Get phenotype information
  temp       <- getPheno.info(pheno.list, snp.list, temp.list=temp.list)
  pheno.id   <- temp$pheno.id
  pheno.uid  <- temp$pheno.uid
  pdata.flag <- temp$pdata.flag
  pdata.ids  <- temp$pdata.ids
  nsubj      <- temp$nsubj
  nsubj.u    <- temp$nsubj.u
  controls   <- getListName(temp, "controls")

  if (op$return.pheno) phenoData.list <- temp$phenoData.list
  rm(pheno.list, temp)
  temp <- gc(verbose=FALSE)

  # Update snp.list
  snp.list <- update.snp.list(snp.list)

  # Define a list to call loadData
  tList <- getSnpLoadList(snp.list, temp.list)

  imputeFlag <- snp.list$file.type %in% c(9, 10)
  if (imputeFlag) {
    op$MAF     <- 0
    op$alleles <- 0
    delimiter  <- "\t"
  }

  snpNames   <- getListName(snp.list, "snpNames")
  total.nsnp <- length(snpNames)
  snpsFlag   <- 1 - is.null(snpNames)
  snp.list$snpNames <- NULL
  count      <- 0
  index      <- 0
  rows       <- NULL
  firstFlag  <- 0
  stopFlag   <- 0
  missing    <- NULL
  snps       <- NULL
  sFlag      <- op$snpNames
  inc.snps   <- op$include.snps
  ret.type   <- op$return.type
  temp.matrix <- NULL
  save.subs  <- NULL
  ret.order  <- NULL
  temp.snp   <- NULL
  MAF.flag   <- op$MAF
  MAF        <- NULL
  allFlag    <- op$alleles
  alleles    <- NULL
  if (!is.null(snp.list$heter.codes)) allFlag <- 0
  gfreqFlag  <- op$genoFreqs
  nmissFlag  <- op$n.miss
  genoFreqs  <- NULL
  n.miss     <- NULL
  temp2      <- NULL
  imputed    <- NULL

  # If the snp names were specified, then set sFlag to 1 so that
  #  the snp names will be kept
  if (snpsFlag) sFlag <- 1
  if (ret.type == 2) sFlag <- 1

  # Function to convert character to numeric
  if (mode != 3) {
    convert.fun <- as.integer
  } else {
    convert.fun <- as.numeric
  }

  out.miss <- snp.list$out.miss
  if (snp.list$file.type == 4) {
    out.sep <- "|"
  } else {
    out.sep <- snp.list$out.delimiter
  }

  # Loop over each file and combine the objects
  for (file in snp.list$file) {
    index <- index + 1

    # Get the input file name
    temp <- paste(snp.list$dir, file, sep="")

    # Update the list
    tList$start.row <- snp.list$start.vec[index]
    tList$stop.row  <- snp.list$stop.vec[index]
    tList$snpNames  <- snpNames

    # Get the snp data
    snpData <- loadData(temp, tList)
    if (imputeFlag) {
      MAF     <- c(MAF, snpData$MAF)
      alleles <- c(alleles, snpData$alleles)
      imputed <- c(imputed, snpData$imputed)
      snpData <- snpData$snpData
    }
    nsnps <- length(snpData)
    if (!firstFlag) {
      if (nsnps <= 1) next
    } else {
      if (!nsnps) next
    }

    # Get the subjects and check for errors
    temp     <- getSubjIds(snpData[1], snpData[2], delimiter, pheno.uid, 
                           controls=controls)
    subjects <- temp$subjects
    subj.ids <- temp$subj.ids
    subjFlag <- temp$subjFlag
    total.nsubjects <- temp$total.nsubjects
    control.ids     <- getListName(temp, "control.ids")

    # Check for controls
    if (!any(control.ids)) {
      #print("##########################################################")
      #print("NOTE: No controls found to determine the minor allele")
      #print("##########################################################")
      control.ids <- NULL
    }

    # Define a matrix for return type = 2
    if (ret.type == 2) temp.matrix <- matrix(data=NA, nrow=nsnps-1, ncol=nsubj)

    # Save the subjects to preserve the order
    if (!firstFlag) {
    
      if (pdata.flag) {
        # The order will be as in the phenotype data
        subj.order  <- pheno.id
        subj.order2 <- getOrder(pheno.id, subjects)

        # The new subject ids will be written out
        subjects    <- pdata.ids
      } else {
        if (op$orderByPheno) {
          subj.order <- pheno.id
        } else {
          subj.order <- subjects
        }
        subj.order2  <- getOrder(subj.order, subjects) 
        subjects     <- subjects[subj.order2]
      }
 
      if (ret.type == 2) colnames(temp.matrix) <- subjects
      if (op$subjIds) save.subs <- subjects
        
      rm(pdata.ids)
    } else {
      subj.order2 <- getOrder(subj.order, subjects)
    }

    # Set the first row
    snpData[1] <- paste(subjects, sep="", collapse=out.sep)

    # For missing vector
    if (missFlag) temp.miss <- rep(FALSE, times=nsnps-1)

    # snp names vector
    if (sFlag) temp.snp <- character(nsnps-1)

    # MAF
    if (MAF.flag) {
      temp.MAF <- numeric(nsnps-1)
    } else {
      temp.MAF <- NULL
    }

    # alleles
    if (allFlag) {
      temp.all <- character(nsnps-1)
    } else {
      temp.all <- NULL
    }

    # Geno freqs
    if (gfreqFlag) {
      temp.gfreq <- character(nsnps-1)
    } else {
      temp.gfreq <- NULL
    }

    # n.miss
    if (gfreqFlag) {
      temp.nmiss <- numeric(nsnps-1)
    } else {
      temp.nmiss <- NULL
    }

    for (i in 2:nsnps) {

      temp <- getVecFromStr(snpData[i], delimiter=delimiter)

      # Remove the snp name
      snp.name <- temp[1]
      temp     <- temp[-1]
      if (sFlag) temp.snp[i-1] <- snp.name

      # Check for error
      lentemp <- length(temp)
      if (lentemp != total.nsubjects) {
        lenstr <- nchar(snpData[i])
        if ((lentemp == total.nsubjects - 1) && ("" %in% in.miss) && 
            (substr(snpData[i], lenstr, lenstr) == delimiter)) {
          # Add a missing genotype
          temp <- c(temp, "")
        } else {
          stop(paste("ERROR with SNP", snp.name))
        }
      }

      # Use the 0, 1, 2 codes
      if (!imputeFlag) {
        temp2 <- recode.geno(temp, in.miss=in.miss, out.miss=out.miss,
              out.genotypes=codes, heter.codes=heter.codes, subset=control.ids)

        # Vector of recoded genotypes
        temp <- temp2$vec
      }

      # Get the MAF
      if (MAF.flag) temp.MAF[i-1] <- getMAF(temp, sub.vec=control.ids, controls=TRUE) 

      if (allFlag) temp.all[i-1] <- temp2$alleles

      # Geno Frequencies
      if (gfreqFlag) temp.gfreq[i-1] <- paste(table(temp, exclude=NA), collapse="|", sep="")

      # n.miss
      if (nmissFlag) temp.nmiss[i-1] <- sum(is.na(temp))

      # Get the correct subjects
      if (subjFlag) temp <- temp[subj.ids]
  
      # Get the correct order
      temp <- temp[subj.order2]

      # Determine if the snp has missing values
      if (missFlag) {
        if (any(is.na(temp))) temp.miss[i-1] <- TRUE
      }

      if (ret.type == 2) {
        temp.matrix[i-1, ] <- convert.fun(temp)
      } else {
        snpData[i] <- paste(temp, sep="", collapse=out.sep)
        if (inc.snps) snpData[i] <- paste(snp.name, out.sep, snpData[i], sep="")
      }
    } # END: for (i in 2:nsnps)

    # Get the rows by searching for the snp names
    if (snpsFlag) {
      snpNames <- update.snpNames(snpNames, temp.snp)
      count    <- count + length(snpData) - 1

      if (!length(snpNames)) stopFlag <- 1
    } # END: if (snpsFlag)

    # Update
    if (ret.type == 2) {
      rownames(temp.matrix) <- temp.snp
      data.obj <- convert.fun(rbind(data.obj, convert.fun(temp.matrix)))
    } else {
      if (firstFlag) snpData <- snpData[-1]
      data.obj  <- c(data.obj, snpData)
    }
    if (missFlag) missing <- c(missing, temp.miss)
    if (sFlag) snps <- c(snps, temp.snp)
    if (MAF.flag) MAF <- c(MAF, temp.MAF)
    if (allFlag) alleles <- c(alleles, temp.all)
    if (gfreqFlag) genoFreqs <- c(genoFreqs, temp.gfreq)
    if (nmissFlag) n.miss <- c(n.miss, temp.nmiss)
    firstFlag <- 1

    if (stopFlag) break
  } # END: for (file in snp.list$file) 

  # Check for error
  if (snpsFlag) {
    if ((total.nsnp != count) && (length(snpNames))) {
      print("Some SNPs were not found")
      if (op$stopOnError) {
        print(snpNames)
        stop("ERROR: The above SNPs were not found")
      }
    }
  }

  rm(subjects, snpNames, tList, subj.ids, pheno.id, pheno.uid, 
    temp.matrix, subj.order, total.nsubjects, temp.MAF, temp.all, temp2,
    temp.gfreq, temp.nmiss)
  temp <- gc(verbose=FALSE)

  if (!op$include.row1) data.obj <- data.obj[-1]

  list(data=data.obj, missing=missing, snpNames=snps, nsubjects=nsubj,
       subjects=save.subs, order=ret.order, phenoData.list=phenoData.list,
       MAF=MAF, alleles=alleles, n.miss=n.miss, genoFreqs=genoFreqs,
       imputed=imputed)

} # END: loadData.type1

# Function to return phenotype data and file id
loadData.stream <- function(snp.list, pheno.list, temp.list, op=NULL) {

  snp.list$stream <- 1
  op <- default.list(op, c("file.index", "useControls"), 
                     list(1, 1), error=c(0, 0))
  index <- op$file.index

  if (index == 1) {
    # Check the lists. pheno.list is checked in getPheno.info
    snp.list  <- check.snp.list(snp.list)
    temp.list <- check.temp.list(temp.list)

    # Update snp.list for snpNames
    snp.list <- update.snp.list(snp.list)

    # Get phenotype information
    pinfo <- getPheno.info(pheno.list, snp.list, temp.list=temp.list)
  } else {
    pinfo <- NULL
  } 
 
  ccvar <- getListName(pheno.list, "cc.var")
  if (is.null(ccvar)) {
    #print("##########################################################")
    #print("WARNING in loadData.stream: pheno.list$cc.var not specified")
    #print("##########################################################")
    op$useControls <- 0
  }
  if (op$useControls == 0) pinfo$controls <- NULL

  # Define a list to call loadData
  tList <- getSnpLoadList(snp.list, temp.list, op=op)

  temp <- paste(snp.list$dir, snp.list$file[index], sep="")

  # Open the data and get the subject ids
  fid <- loadData(temp, tList)

  list(phenoData.list=pinfo$phenoData.list, pheno.id=pinfo$pheno.id,
      pheno.uid=pinfo$pheno.uid, fid=fid, controls=pinfo$controls)

} # END: loadData.stream

# Function to get the data or the file id
getData.1 <- function(snp.list, pheno.list, temp.list, op=NULL) {

  # Check the lists
  snp.list   <- check.snp.list(snp.list)
  pheno.list <- check.pheno.list(pheno.list)

  # Check the delimiters
  temp <- snp.list
  temp$file <- temp$file[1]
  if (!checkDelimiter(temp)) {
    stop("ERROR: Check the delimiter in snp.list")
  }
  rm(temp)
  if (!checkDelimiter(pheno.list)) {
    stop("ERROR: Check the delimiter in pheno.list")
  }

  if (snp.list$stream) {
    ret <- loadData.stream(snp.list, pheno.list, temp.list, op=op)
  } else {
    ret <- loadData.type1(snp.list, pheno.list, temp.list, op=op)
  }

  ret

} # END: getData.1

# Function to get the next observation 
scanNextObs <- function(fid, fid.row, get.row, snpFlag=0, sep="|",
                      snpNames=NULL) {
 
  if (!snpFlag) {
    row <- scan(file=fid, what="character", nlines=1, sep=sep,
                      skip=get.row-fid.row, quiet=TRUE)
  } else {
    # Loop and check each snp name
    cont  <- 1
    while (cont) {
      row <- scan(file=fid, what="character", nlines=1, sep=sep, quiet=TRUE)
      if (!length(row)) break
      if (row[1] %in% snpNames) break 
    }
  }

  row

} # END: scanNextObs

# Function to get the next observation and close the file if no more
#   reading is to be done
getNextObs <- function(i, snpFid, snpFlag, snpNames, tempfile, delete,
                  start.stream, stop.stream, delimiter) {

  # Check the snpNames vector
  if ( ((snpFlag) && (!length(snpNames))) || (i > stop.stream) ) {
    closeFile(snpFid, file=tempfile, delete=delete)
    return(NULL)
  }

  # Read in the next snp
  snp <- scanNextObs(snpFid, start.stream, start.stream, 
           sep=delimiter, snpFlag=snpFlag, snpNames=snpNames)
  if (!length(snp)) {
    closeFile(snpFid, file=tempfile, delete=delete)
    return(NULL)
  }
     
  snp

} # END: getNextObs

# Function to order and subset a SNP
orderSNP <- function(snp, snp.name, subj.order2, total.nsubjects, 
                     in.miss, heter.codes, subjFlag=0, subj.ids=NULL, 
                     out.genotypes=0:2, control.ids=NULL) {

  # Check for error
  if (length(snp) != total.nsubjects) stop(paste("ERROR with SNP", snp.name))
 
  # Use the 0, 1, 2 codes
  temp <- recode.geno(snp, in.miss=in.miss, out.miss=NA,
                 out.genotypes=out.genotypes, heter.codes=heter.codes,
                 subset=control.ids)
  snp <- temp$vec
  alleles <- temp$alleles

  # Compute the MAF
  maf <- getMAF(snp, sub.vec=control.ids, controls=TRUE)

  # Get the correct subjects
  if (subjFlag) snp <- snp[subj.ids]
  
  # Get the correct order
  snp <- snp[subj.order2]

  list(SNP=snp, MAF=maf, alleles=alleles)

} # END: orderSNP

# Function to call for stream input after getData.1 was called
setUp.stream <- function(obj, snp.list, tempfile, delete, controls=NULL) {
  
  # obj        Return object from getData.1
  
  pheno.uid <- obj$pheno.uid
  pheno.id  <- obj$pheno.id
  snpFid    <- obj$fid
  start.vec <- snp.list[["start.vec", exact=TRUE]]
  stop.vec  <- snp.list[["stop.vec", exact=TRUE]]
  start     <- max(2, start.vec[1])
  # Get stop in terms of i = 1, 2, ....
  stop      <- stop.vec[1] - start + 1
  if (stop < 0) stop <- Inf
  delimiter <- getDelimiter(snp.list)
  snames    <- snp.list[["snpNames", exact=TRUE]]
  snpFlag   <- 1 - is.null(snames)

  # Read in the subjects and first snp
  temp <- scanNextObs(snpFid, 1, 1, sep=delimiter)
  snp  <- scanNextObs(snpFid, 2, start, sep=delimiter, 
                     snpFlag=snpFlag, snpNames=snames)
  if (!length(snp)) {
     closeFile(snpFid, file=tempfile, delete=delete)
     stop("ERROR: Check snp.list options")
  }
 
  # Get the subjects id vector and order
  temp            <- getSubjIds(temp, snp, delimiter, pheno.uid, controls=controls)
  subjects        <- temp$subjects
  subjFlag        <- temp$subjFlag
  subj.ids        <- temp$subj.ids
  total.nsubjects <- temp$total.nsubjects
  control.ids     <- getListName(temp, "control.ids")
  if (!any(control.ids)) {
    #print("##########################################################")
    #print("NOTE: No controls found to determine the minor allele")
    #print("##########################################################")
    control.ids <- NULL
  }
  subj.order2     <- getOrder(pheno.id, subjects)

  list(subjFlag=subjFlag, subj.ids=subj.ids, start.stream=start,
       total.nsubjects=total.nsubjects, subj.order2=subj.order2,
       snp=snp, stop.stream=stop, control.ids=control.ids)

} # END: setUp.stream

# Function to get info from the phenotype data
getPheno.info <- function(pheno.list, snp.list, temp.list=NULL) {

  id.var <- pheno.list$id.var

  # Load the phenotype data
  temp <- list(file.type=pheno.list$file.type, header=pheno.list$header,
               delimiter=pheno.list$delimiter, id.var=id.var,
               remove.miss=0, in.miss=pheno.list$in.miss)
  data <- loadData(pheno.list$file, temp)

  # Get all of the control ids
  ccvar <- pheno.list[["cc.var", exact=TRUE]]
  if (!is.null(ccvar)) {
    # Check that it is 0-1
    temp <- (data[, ccvar] %in% c(0, 1))
    temp[is.na(temp)] <- TRUE
    if (!all(temp)) stop("ERROR in getPheno.info: pheno.list$cc.var is not coded as 0-1")

    temp     <- (data[, ccvar] == 0)
    temp[is.na(temp)] <- FALSE
    controls <- unique(makeVector(data[temp, id.var]))
  } else {
    controls <- NULL
  }

  # Update pheno.list
  pheno.list$data        <- data
  pheno.list$is.the.data <- 1

  # Get the intersecting subject ids with the genotype data
  ids <- intersectSubIds(snp.list, pheno.list, temp.list=temp.list)

  # Keep these ids
  temp <- data[, id.var] %in% ids
  temp[is.na(temp)] <- FALSE
  data <- removeOrKeepRows(data, temp)
  pheno.list$data <- data
  rm(data, ids)
  gc()

  # Call getPhenoData for the other pheno.list options to get the 
  # data to be used in the analysis
  getPheno.list      <- getPhenoData(pheno.list, temp.list=temp.list)
  phenoData          <- getPheno.list$data
  pheno.list$data    <- NULL
  pdata.flag         <- !is.null(getPheno.list$orig.id)
  pheno.list$orig.id <- getPheno.list$orig.id

  # pdata.flag is the flag for non-unique subject ids
  if (!pdata.flag) {
    # Get the subject ids
    pheno.id  <- makeVector(phenoData[, pheno.list$id.var])
    pdata.ids <- pheno.id
  } else {
    pheno.id  <- as.character(phenoData[, pheno.list$orig.id])
    pdata.ids <- makeVector(phenoData[, pheno.list$id.var])
  }

  # Get the total number of subjects
  nsubj <- length(pheno.id)

  # Get the unique ids
  pheno.uid <- unique(pheno.id)
   
  # Get the number of unique subjects
  nsubj.u <- length(pheno.uid)

  # Return a list of info
  list(pheno.id=pheno.id, pheno.uid=pheno.uid, pdata.ids=pdata.ids,
       pdata.flag=pdata.flag, nsubj=nsubj, nsubj.u=nsubj.u,
       phenoData.list=getPheno.list, controls=controls)

} # END: getPheno.info

# Function to check for errors is the subject ids
getSubjIds <- function(snpData1, snpData2, delimiter, pheno.uid, controls=NULL) {

  # snpData1      Row 1      
  # snpData2      Row 2
  # delimiter     delimiter
  # pheno.uid     Unique (original) phenotype ids
  # controls      Unique (original) control ids

  # Get the subject ids from the snp data
  subs <- getSubject.vec(snpData1, snpData2, delimiter)
  total.nsubjects <- length(subs)

  # Get the number of unique original subject ids
  nsubj.u <- length(pheno.uid)

  # Check for an error
  if (sum(pheno.uid %in% subs) != nsubj.u) {
    temp <- !(pheno.uid %in% subs)
    temp <- pheno.uid[temp]
    #bad.subject.ids.global <<- temp
    print(temp)
    #print("bad.subject.ids.global")
    stop("The above subject ids were not found in the genotype data")
  }

  # Determine if any subjects should be removed
  subj.ids <- subs %in% pheno.uid

  if (sum(subj.ids) != nsubj.u) {
    print("The subject ids in the genotype data may not all be unique")
    stop("ERROR with subject ids")
  }
  
  # Get the logical vector for controls. Recall: MAF is determined by the 
  #  entire sample of subjects.
  if (!is.null(controls)) {
    control.ids <- subs %in% controls
  } else {
    control.ids <- NULL
  }

  if (total.nsubjects == nsubj.u) {
    # All subjects are here, so subj.ids is no longer needed
    subj.ids <- NULL
    subjFlag <- 0
  } else {
    subjFlag <- 1
    subs     <- subs[subj.ids]
  }

  # Return list
  list(subjects=subs, subjFlag=subjFlag, subj.ids=subj.ids,
       total.nsubjects=total.nsubjects, control.ids=control.ids)

} # END: getSubjIds

# Function to return the vector of subject ids from a snp file
getSubject.vec <- function(snpData1, snpData2, delimiter) {

  # snpData1    Row 1 (header) 
  # snpData2    Row 2
  # delimiter

  # Get the subject ids from the snp data
  if (length(snpData1) == 1) {
    subs <- getVecFromStr(snpData1, delimiter=delimiter)
  } else {
    subs <- snpData1
  }

  # Get the number of fields to remove
  n.omit <- getRow1.omitN(subs, snpData2, delimiter=delimiter) 
  if (n.omit > 0) subs <- subs[-(1:n.omit)] 

  subs

} # END: getSubject.vec

# Function to determine the number of fields to remove in row1 of the snp
#  data 
getRow1.omitN <- function(row1, row2, delimiter="|") {

  # row1  
  # row2  
  # delimiter

  n1 <- length(row1)
  if (length(row2) == 1) row2 <- getVecFromStr(row2, delimiter=delimiter)
  n2 <- length(row2)

  n.omit <- n1 - n2 + 1
  if (n.omit < 0) stop("ERROR: in the snp files")
  if (n.omit > 1) warning("Possible error in the SNP data")

  n.omit
 
} # END: getRow1.omitN

# Function to return the delimiter used in the file read in
getDelimiter <- function(snp.list, output=0) {

  delimiter <- snp.list$delimiter
  if ((output == 1) && (!snp.list$stream)) delimiter <- snp.list$out.delimiter
  if (snp.list$file.type %in% c(4)) {
    delimiter <- "|"
  } 
  delimiter

} # END: getDelimiter

# Function to get the vector of missing values
get.in.miss <- function(snp.list) {

  ret <- snp.list$in.miss
  if (snp.list$file.type == 4) ret <- c(ret, "-9")
  ret

} # END: get.in.miss

# Function to update the snpNames vector when snp names are specified
update.snpNames <- function(snpNames, temp.snp) {

  # snpNames    Vector of all snp names desired
  # temp.snp    Vector of snp names on current data

  # Remove the snp names from the snpNames vector
  rows     <- as.logical(1 - (snpNames %in% temp.snp))
  snpNames <- snpNames[rows]
  
  snpNames

} # END: update.snpNames

# Function to return a list of parameter options for loading the 
#  snp data.
getSnpLoadList <- function(snp.list, temp.list, op=NULL) {

  op <- default.list(op, c("file.index"), list(1))

  ret <- list(file.type=snp.list$file.type, 
         delimiter=snp.list$delimiter, read.n=snp.list$read.n,
         sas.list=snp.list$sas.list, transpose=1, include.row1=1,
         id.var=snp.list$id.var, 
         start.row=snp.list$start.vec[op$file.index],
         stop.row=snp.list$stop.vec[op$file.index],
         snpNames.keep=snp.list$snpNames.keep)
  ret$snpNames <- getListName(snp.list, "snpNames")
  stream <- snp.list$stream
  if (stream) {
    ret$stream  <- 1
    ret$outfile <- op$outfile
  }
  if (ret$file.type == 4) ret$sas.list$temp.list <- temp.list
  if (ret$file.type %in% c(9, 10)) {
    ret <- snp.list
    ret$delimiter <- "\t"
    ret$out.delimiter <- "\t"
  }
  
  ret

} # END: getSnpLoadList

# Function to get the common subject ids.
# Returns the updated pheno.list
intersectSubIds <- function(snp.list, pheno.list, temp.list=NULL) {

  snp.list$snpNames      <- NULL
  snp.list$snpNames.list <- NULL
  snp.list   <- check.snp.list(snp.list)
  pheno.list <- check.pheno.list(pheno.list)
  temp.list  <- check.temp.list(temp.list)
  delimiter  <- getDelimiter(snp.list)
  snp.list   <- update.snp.list(snp.list)
  stream     <- getListName(snp.list, "stream")

  # Get the unique files
  files <- unique(snp.list$file)

  # Define a list to call loadData
  tList           <- getSnpLoadList(snp.list, temp.list)
  tList$start.row <- 1
  tList$stop.row  <- 2
  isub            <- NULL
  index           <- 1

  # Loop over each file and combine the objects
  for (file in files) {

    # Get the data
    temp <- loadData(paste(snp.list$dir, file, sep=""), tList)
    if (snp.list$file.type %in% c(9, 10)) {
      temp <- temp$snpData
      delimiter <- "\t"
    }

    # Get the first 2 rows
    if (!stream) {
      row1 <- temp[1]
      row2 <- temp[2]
    } else {
      row1 <- scanNextObs(temp, 1, 1, sep=delimiter, snpFlag=0, 
                          snpNames=NULL)
      row2 <- scanNextObs(temp, 1, 1, sep=delimiter, snpFlag=0, 
                          snpNames=NULL)
      close(temp)
    }

    # Get the subject ids
    temp <- getSubject.vec(row1, row2, delimiter)
    if (index == 1) {
      isub <- temp
    } else {
      isub <- intersect(isub, temp)
    }
  }

  if (!length(isub)) {
    stop("ERROR: No intersecting subject ids in the genotype files")
  }

  # Get the phenotype data
  dflag <- !is.null(getListName(pheno.list, "is.the.data"))
  if (dflag) {
    temp <- unique(makeVector(pheno.list$data[, pheno.list$id.var]))
  } else {
    temp <- NULL # ADD code here
  }

  isub <- intersect(isub, temp)
  if (!length(isub)) {
    stop("ERROR: No intersecting subject ids in the genotype and phenotype files")
  }

  if (length(temp) != length(isub)) {
    temp <- temp[!(temp %in% isub)]
    print(temp)
    #bad.subject.ids.global <<- temp
    print("The above subject ids were not found in the genotype data")
    warning("Some subject ids were not found in the genotype data")
  }

  isub

} # END: intersectSubIds
