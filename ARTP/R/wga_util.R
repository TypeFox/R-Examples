# History  Mar 28 2008 Add option to createDummy
#                      Add getOrder function
#          Mar 31 2008 Change snpNames option. snpNames are added
#                      in check.snp.list
#          Apr 07 2008 Add getNames function.
#                      Check for file seperator in check list functions.
#          Apr 11 2008 Add function removeOrKeepCols
#          Apr 15 2008 Change to sas.list
#          Apr 16 2008 Add removeOrKeepRows function
#          Apr 17 2008 Change in default.list
#          Apr 25 2008 Add unfactor function
#          Apr 30 2008 Add getVarNames function
#          May 01 2008 Add factorVars function
#                      Allow logical vector in removeOrKeepRows
#          May 05 2008 Change getSnpNames to getIdsFromFile
#          May 06 2008 Change snp.col, chrm.col, and loc.col to
#                      snp.var, chrm.var, and loc.var
#          May 09 2008 Add getInheritanceVec function
#                      Add addInterVars function
#          May 13 2008 Add option to createDummy to keep the original
#                      factor variable in the returned data frame.
#          May 31 2008 Add function matchNames
#          Jun 07 2008 Unfactor factors in changeStrata
#          Jun 10 2008 Change inheritance to genetic.model
#          Jun 11 2008 Add checkList option to default.list function
#          Jun 13 2008 Add getCommandArg function
#          Jun 14 2008 Remove nvars in getVarNames
#          Jun 20 2008 Add function to define lists for genfile
#                      Make snp.list$in.miss and snp.list$heter.codes
#                      to be character vectors in sheck.snp.list
#                      Check if input list is NULL in default.list
#                      Fix bug in recode.geno
#          Jun 21 2008 Change createDummy to allow for input matrices
#                      Return a list with added variable names
#          Jun 25 2008 Make createDummy more efficient
#                      Add start.n to changeStrata
#                      Make removeMiss more general
#                      Add sort2D function
#          Jun 27 2008 Change createDummy function
#          Jun 27 2008 Use getColumns function
#          Jun 27 2008 Add update.snp.list function
#                      Add getTempfile
#          Jun 28 2008 Add function extractByStr
#          Jun 30 2008 Add function closeFile
#          Jul 02 2008 Add code for type 5 files in check.snp.list
#                      Allow for GLU formats
#                      Add makeVector function
#          Jul 06 2008 Add initDataFrame function
#          Jul 08 2008 Add functions for subsetting data
#          Jul 11 2008 Add renameVar function
#          Jul 15 2008 Make default delimiter in snp.list a tab
#          Jul 17 2008 Change addInterVars
#                      Add getVarNames.int
#          Jul 24 2008 Add getSNP.type2
#          Jul 30 2008 Change to checkForSep, check.snp.list, 
#                       check.locusMap.list
#          Aug 07 2008 Add checkForConstantVar function
#          Aug 11 2008 Add unfactor.all function
#          Aug 15 2008 Redo getSNPdata
#          Aug 18 2008 Add strataMatrix
#          Aug 19 2008 Add rep.rows, rep.cols, rep.mat
#          Oct 10 2008 Add recode.data function
#          Oct 11 2008 Fix bug in subsetData.var for vector value
#          Oct 21 2008 Fix bug in addInterVars
#          Oct 22 2008 Extend recode.data function
#                      Add mergeData function
#          Oct 29 2008 Make changeLevels.var more general
#          Oct 31 2008 Allow for user defined genotypes in recode.data
#          Nov 04 2008 Extend removeOrKeepRows for character rows
#          Nov 05 2008 Extend subsetData.var to return rows
#          Nov 06 2008 Change the functionality subsetData.list
#                      Extend subsetData.var for character vector values
#          Nov 12 2008 Add callOS function
#          Nov 26 2008 Add GetColsFromCharVec function
#          Nov 26 2008 Generalize removeMiss function
#          Dec 08 2008 Add orderGenotypes function
#          Dec 10 2008 Keep rownames in removeOrKeepCols/removeOrKeepRows
#                      Add changeStr.names function
#          Dec 11 2008 Add applyFormulas function
#                      Add removeMiss.vars function
#          Jan 08 2009 Add changeAlleles function
#          Jan 21 2009 Add addColumn function
#          Jan 30 2009 Update changeAlleles function.
#                      Add writeTable function
#          Jan 31 2009 Update subsetData.var for missing values
#          Feb 02 2009 Add getPartition function
#          Feb 04 2009 Add crossTab function
#          Feb 05 2009 Change in recode.geno to allow for determining the major
#                      and minor alleles from a subset. 
#          Feb 07 2009 Add mergePhenoGeno function
#          Feb 23 2009 Fix bug in subsetData.var with NAs
#          Feb 26 2009 Extend crossTab function
#          Mar 05 2009 Make addColumn more efficient
#          Mar 16 2009 Update crossTab for 1 variable
#                      Add parse.vec function
#          Mar 19 2009 Add error checks in subsetData.list
#          Mar 20 2009 Change to subsetData.list to return the rows
#          Mar 23 2009 Update sort2D for 1 NA row
#          Mar 24 2009 Add replaceStr.var
#          Mar 25 2009 Add warnings to addColumn function
#          Mar 27 2009 Add option returnRows to subsetData.list
#          Apr 01 2009 Add option for NAs in subsetData
#          Apr 07 2009 Fix bug in recode.geno
#          Apr 08 2009 Let file be a data frame in addColumn
#          Apr 14 2009 Generalize subsetData.list
#          May 07 2009 Add option to mergePhenoGeno to include original genotypes
#          May 07 2009 Change replaceStr.var function
#                      Combine addColumn and addCol2
#                      Remove addCol2 function
#          May 12 2009 Add function orderVars 
#          May 13 2009 Add option to addColumn to check for leading zeros
#                      in the id variables.
#          May 29 2009 Update removeOrKeepCols to print out cols not found
#          Jun 01 2009 Update addColumn
#          Jun 05 2009 Add getTopHits function
#          Jun 05 2009 Add option to addColumn for not replacing variables.
#                      Add function mergeTopHits
#          Jun 08 2009 Add getRank.file function
#          Jun 09 2009 Change in subsetData.var for a numeric variable
#          Jun 15 2009 Add option to extractByStr for efficiency
#          Jun 25 2009 Force input argument data to be a data frame in addColumn
#                      Create function addRow
#          Jul 02 2009 Return a matrix in addColumn if the input object is a matrix
#          Jul 02 2009 Add option to sort2D
#          Jul 10 2009 Fix bug in getVarNames with NULL input object
#          Jul 21 2009 Add option to removeOrKeepCols to not throw error
#                      when columns are not found
#                      Move out older code to wga_util2.R
#                      Add function to check the delimiter in a file
#          Jul 22 2009 Add flag to snp.list, pheno.list, etc to check if
#                      the check.xxx.list was called.
#                      Add name mergeDirFile to snp.list
#          Jul 24 2009 Add writeVec function
#          Aug 03 2009 Changes in check.pheno.list
#          Aug 04 2009 Add removeWhiteSpace function
#                      Move less used code to wga_util2.R
#          Aug 10 2009 Check snpNames.list in check.snp.list
#          Aug 10 2009 Add options to snp.list in check.snp.list
#                      Add option in extractByStr to keep or remove
#                      matched values
#          Aug 18 2009 Add debug.time function
#          Aug 19 2009 Update check.vec to return flag
#          Aug 27 2009 Add exclude option in crossTab
#          Sep 24 2009 Add function to compute all interactions
#                      between 2 sets of variables in a data frame
#          Oct 01 2009 Change var to "var" in checkForConstantVar function
#          Oct 16 2009 Set up changes for efficiency:
#                      check.snp.list
#                      check.pheno.list
#          Oct 16 2009 Return major/minor alleles in recode.geno
#          Oct 16 2009 Change in closeFile to check that the file exists
#                      before deleting it.
#          Oct 20 2009 Add fun option to unfactor.all
#          Oct 23 2009 Add function checkTryError
#          Oct 27 2009 Change else statement in getInheritanceVec
#          Oct 28 2009 Add option removeVars from checkForConstantVar
#          Oct 28 2009 Fix bug in checkForConstantVar
#          Nov 03 2009 Add column names in getColsFromCharVec
#          Nov 10 2009 Update genfile.list for functions
#          Nov 10 2009 Update crossTab so that snp.list can be NULL
#          Nov 12 2009 Fix bug in recode.geno when out.genotypes is NULL
#          Nov 12 2009 Check vars in addColumn
#          Nov 20 2009 NO CHANGE
#          Dec 10 2009 Fix bug in determining the major/minor alleles in recode.geno
#          Dec 17 2009 Add option to getVarNames.int for seperating 
#                       interaction terms      
#          Dec 22 2009 Initialize a1, a2 to NULL in recode.geno       
#          Dec 29 2009 Use getFileDelim in check.xxx.list 
#          Dec 30 2009 Add code for getting variables from formulas    
#                      Add miss option to removeMiss.vars      
#                      Set default value for in.miss in check.pheno.list 
#          Dec 31 2009 Add option to callOS
#          Jan 04 2010 Update applyFormulas
#          Jan 13 2010 Add function check.file.list
#          Jan 15 2010 Add check.subsetData.list, update check.file.list.
#                      Add getSubsetDataVars, check.subsetData.list
#          Jan 19 2010 Update check.pheno.list with check.file.list
#          Feb 01 2010 Set extended = FALSE in replaceStr.var
#          Mar 04 2010 Add normVarNames function
#          Mar 05 2010 Update getIdsFromFile to call loadData.table
#          Mar 14 2010 Check vector length in recode.geno
#          Mar 18 2010 Fix bug in genfile.list with an empty list
#          Apr 13 2010 Update mergePhenoGeno for imputed data
#          Apr 15 2010 Add allele option to mergePhenoGeno
#          Apr 27 2010 Add new pheno.list options for duplicated ids
#          May 04 2010 No changes
#          Jun 14 2010 Add return option to getColsFromCharVec
#          Jul 02 2010 Add option to return snp names in mergePhenoGeno
#                      Update writeTable function
#          Jul 22 2010 Update addColumn
#          Jan 04 2011 Update check.file.list with getFileHeader function
#          Jan 06 2011 Update extractByStr to allow for exact matching
#          Jan 12 2011 Remove extended option in grep, gsub
#          Jan 25 2011 Add changeGenotypes function
#                      Add getSnpNames function
#          Jan 31 2011 Allow snp.list and pheno.list to be file names instead of lists
#                      Remove all calls to getListName
#          Jul 12 2011 Add function to return a unique variable name
#          Dec 22 2011 Allow for genetic.model = 4 (heterozygous) in check.snp.list
#          Jan 13 2012 No change
#          Feb 09 2012 Add code in recode.geno to check heterozygous codes when more than 3 genotypes

# Function to pull apart each string from the data object
getVecFromStr <- function(string, delimiter="|") {

  # string       String to break apart. No default
  # delimiter    Delimiter used in string. The default is "|". 

  strsplit(string, delimiter, fixed=TRUE)[[1]]

} # END: getVecFromStr

# Function to recode a vector of genotypes and missing values
recode.geno <- function(vec, in.miss=c("  "), out.miss=NA,
               out.genotypes=c(0,1,2), heter.codes=NULL, subset=NULL) {

  # vec           Input vector 
  #               No default.
  # in.miss       Vector of categories which denote missing values. 
  #               These categories will be changed to out.miss,
  #               if out.miss is not NULL.
  #               The default is "  "  (2 spaces)
  # out.miss      New category for missing values. Use NULL if you want to
  #               keep the original missing categories.
  #               The default is NA
  # out.genotypes Vector of length 3 containing the new genotypes,
  #               the first element is for the major homozygous
  #               genotype and the second element is for the 
  #               heterozygous genotype. 
  #               Use NULL for no recoding of the genotypes
  # heter.codes   Vector of codes used for the heterozygous genotype.
  #               If NULL, then it is assumed that the heterozygous 
  #               genotype is of the form "AB", "Aa", "CT", etc (a 2-character
  #               string with different characters... case sensitive!!!)
  #               The default is NULL.
  # subset        NULL or a logical vector of length length(vec) to be used
  #               in determining the major and minor homozygous genotypes.
  #               The default is NULL so that all (non-missing) elements
  #               of vec will be used. 

  # Get the rows that contain missing values
  # Do not change these values now, because we need to know the other codes.
  tempMiss <- vec %in% in.miss
  if (!is.null(out.miss)) {
    mFlag    <- 1
  } else {
    mFlag    <- 0
  }

  # Determine if a subset is to be used
  subFlag <- !is.null(subset)
  if (subFlag) {
    if (!is.logical(subset)) stop("ERROR in recode.geno: subset is not a logical vector")
    if (length(subset) != length(vec)) stop("ERROR in recode.geno: length(subset) != length(vec)")
    if (!any(subset)) subFlag <- 0
    vec0 <- vec
  }

  # Initialize
  index   <- Inf
  subSum  <- Inf
  alleles <- "  "
  hflag   <- 1
  a1      <- NULL
  a2      <- NULL

  if (!is.null(out.genotypes)) {

    # Watch out for heterozygous codes. If snp vector contains more than 1 distinct
    # heterozygous code, then make them consistent
    if (!is.null(heter.codes)) {
      temp <- vec %in% heter.codes
      if (any(temp)) {
        uht  <- unique(vec[temp]) 
        nuht <- length(uht)
        if (nuht > 1) vec[temp] <- uht[1] 
      }
    }

    # Get the frequency counts
    if (subFlag) {
      tab <- sort(table(vec[subset], exclude=in.miss), decreasing=TRUE)
    } else {
      tab <- sort(table(vec, exclude=in.miss), decreasing=TRUE)
    }

    # Remove the missing values
    if (!is.null(in.miss)) {
      genos <- names(tab)
      tab   <- tab[!(genos %in% in.miss)]
    }

    # Check for error
    if (length(tab) > 3) {
      print(tab)
      stop("ERROR: in recode.geno")
    }

    # Get the genotypes
    genos   <- names(tab)
    flag    <- 0
    flag2   <- 0
    minFlag <- 0
    hflag   <- !is.null(heter.codes)  

    # Get the ids for each genotype
    ids    <- list()
    index  <- 1
    if (subFlag) subSum <- sum(tempMiss)
    for (geno in genos) {
      ids[[index]] <- vec == geno
      if (subFlag) subSum <- subSum + sum(ids[[index]])
      index        <- index + 1   
    } 

    # tab is sorted in descending order
    index <- 1
    for (geno in genos) {

      # Determine if this genotype is heterozygous
      if (hflag) {
        heterFlag <- geno %in% heter.codes
      } else {
        a1        <- substr(geno, 1, 1)
        a2        <- substr(geno, 2, 2)
        heterFlag <- a1 != a2
      }

      if (!heterFlag) {
        # Homozygous, but check if major homozygous has been assigned
        if (!flag) {
          vec[ids[[index]]] <- out.genotypes[1]
          flag <- 1
          hom1 <- a1
        } else {
          # Minor homozygous
          vec[ids[[index]]] <- out.genotypes[3]
          minFlag <- 1
          hom2    <- a1
        }
      } else {
        # Heterozygous
        # Check for error 
        if (flag2) {
          print(tab)
          stop("ERROR: in recode.geno, flag2=1")
        }
        vec[ids[[index]]] <- out.genotypes[2]
        flag2 <- 1
        hetg  <- geno
      }
      index <- index + 1

    } # END: for (geno in genos)

  } # END: if (!is.null(out.genotypes))

  # Change missing values
  if (mFlag) vec[tempMiss] <- out.miss

  # Get the major/minor alleles
  if (!hflag) {
    if (flag) {
      # Major
      if (minFlag) {
        # There was major and minor genotypes
        alleles <- paste(hom1, hom2, sep="")
      } else if (flag2) {
        # Use heterozygous genotype
        a2 <- substr(hetg, 2, 2)
        if (hom1 != a2) {
          alleles <- paste(hom1, a2, sep="")
        } else {
          a1 <- substr(hetg, 1, 1)
          alleles <- paste(hom1, a1, sep="")
        }
      } else {
        # Only major homozygous
        alleles <- paste(hom1, hom1, sep="") 
      }
    } else {
      # Only heterozygous genotype
      if (flag2) alleles <- hetg
    }
  }

  # Final error check if a subset was used
  if ((subFlag) && (index < 4) && (subSum < length(vec))) {
    temp <- as.numeric(!flag) + as.numeric(!flag2) + as.numeric(!minFlag)
    if (temp > 1) {
      # Use all elements
      temp <- recode.geno(vec0, in.miss=in.miss, out.miss=out.miss,
               out.genotypes=out.genotypes, heter.codes=heter.codes, subset=NULL)
      return(temp)
    }

    # Get the ids that were not mapped
    temp <- tempMiss
    for (i in 1:length(ids)) temp <- temp | ids[[i]]
      
    if (!minFlag) vec[!temp] <- out.genotypes[3]
    if (!flag2)   vec[!temp] <- out.genotypes[2]
    if (!flag)    vec[!temp] <- out.genotypes[1]
  }

  list(vec=vec, alleles=alleles)

} # END: recode.geno

# Function to assign a default value to an element in a list
default.list <- function(inList, names, default, error=NULL,
                         checkList=NULL) {

  # inList      List
  # names       Vector of names of items in inList
  # default     List of default values to assign if a name is not found
  #             The order of default must be the same as in names.
  # error       Vector of TRUE/FALSE if it is an error not to have the
  #             name in the list. 
  #             The default is NULL
  # checkList   List of valid values for each name.
  #             Use NA to skip a list element.
  #             The default is NULL

  n1 <- length(names)
  n2 <- length(default)
  if (n1 != n2) stop("ERROR: in calling default.list")

  if (is.null(error)) {
    error <- rep(0, times=n1)
  } else if (n1 != length(error)) {
    stop("ERROR: in calling default.list")
  }

  if (!is.null(checkList)) {
    if (n1 != length(checkList)) stop("ERROR: in calling default.list")
    checkFlag <- 1
  } else {
    checkFlag <- 0
  } 

  if (is.null(inList)) inList <- list()

  listNames <- names(inList)
  for (i in 1:n1) {
    if (!(names[i] %in% listNames)) {
      if (!error[i]) {
        inList[[names[i]]] <- default[[i]]
      } else {
        temp <- paste("ERROR: the name ", names[i], " was not found", sep="")
        stop(temp)
      }
    } else if (checkFlag) {
      temp <- checkList[[i]]
      if (!all(is.na(temp))) {
        if (!all(inList[[names[i]]] %in% checkList[[i]])) {
          temp <- paste("ERROR: the name '", names[i], 
                      "' has an invalid value", sep="")
          stop(temp)
        }
      }
    }
  }

  inList

} # END: default.list

# Function to get a name from a list (without partial matching)
getListName <- function(inList, name) {

  if (name %in% names(inList)) {
    return(inList[[name]])
  } else {
    return(NULL)
  }

} # END: getListName

# Function to check snp.list
check.snp.list <- function(snp.list) {

  if (is.null(snp.list)) stop("ERROR: snp.list must be specified")
  if (is.character(snp.list)) snp.list <- list(file=snp.list)

  # Check the names in the list
  snp.list <- default.list(snp.list, 
            c("file", "in.miss", 
        "read.n", "genetic.model", "stream", "recode",
        "alreadyChecked", "out.miss", "out.delimiter", "snpNames.keep"), 
        list("ERROR", "  ", -1, 0, 0, 1, 0, NA, "\t", 1), 
            error=c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0),
         checkList=list(NA, NA, NA, 0:4, 0:1, 0:1, NA, NA, NA, NA))

  if (snp.list$alreadyChecked == 1) return(snp.list)

  # Check file vector
  nfiles <- length(snp.list$file)
  snp.list$dir <- checkForSep(snp.list[["dir", exact=TRUE]])
  temp <- paste(snp.list$dir, snp.list$file, sep="")
  if (check.files(temp)) stop("ERROR: in check.snp.list") 
  snp.list$file <- temp
  snp.list$dir <- ""

  if (is.null(snp.list[["file.type", exact=TRUE]])) {
    snp.list$file.type <- getFileType(snp.list$file[1], default=7)
  }
  if (is.null(snp.list[["delimiter", exact=TRUE]])) {
    snp.list$delimiter <- getFileDelim(snp.list$file[1], type=snp.list$file.type, default="\t")
  }

  # Check start.vec and stop.vec
  snp.list <- default.list(snp.list, c("start.vec", "stop.vec"), 
               list(rep(1, nfiles), rep(-1, nfiles)), error=c(0, 0))
  temp <- snp.list[["snpNames", exact=TRUE]]
  if (is.null(temp) && is.null(snp.list$snpNames.list)) {
    a1 <- check.vec.num(snp.list$start.vec, "snp.list$start.vec", len=nfiles)
    a2 <- check.vec.num(snp.list$stop.vec, "snp.list$stop.vec", len=nfiles)
    if (sum(a1+a2)) stop() 
  } else {
    snp.list$start.vec <- rep(1, times=nfiles)
    snp.list$stop.vec  <- rep(-1, times=nfiles)
  }
  temp <- (snp.list$stop.vec == 1)
  if (any(temp)) snp.list$stop.vec[temp] <- 2
  
  # Check for id variable for type 3 and 4
  if (snp.list$file.type %in% c(3, 4)) {
    #if (is.null(snp.list$id.var)) stop("snp.list$id.var is not specified")
    if (is.null(snp.list$id.var)) snp.list$id.var <- 1
  }

  # Check for sas.list
  if (snp.list$file.type == 4) {
    if (is.null(snp.list$sas.list)) stop("snp.list$sas.list is not specified")
  }

  # Check for zip file with type 5
  if (snp.list$file.type == 5) {
    if (is.null(snp.list$zipFile)) stop("snp.list$zipFile is not specified")
  }

  # Check in.miss
  vec  <- as.character(snp.list$in.miss)
  temp <- is.na(vec)
  if (any(temp)) vec[temp] <- "NA"
  snp.list$in.miss <- vec

  # Check heter.codes
  vec <- snp.list$heter.codes
  if (!is.null(vec)) {
    vec <- as.character(vec)
    temp <- is.na(vec)
    if (any(temp)) vec[temp] <- "NA"
    snp.list$heter.codes <- vec
  }

  # Check file.type
  temp <- snp.list$file.type
  if (is.numeric(temp)) {
    if (!(temp %in% 1:10)) {
      temp <- paste("ERROR:", temp, "is not a valid value for snp.list$file.type")
      stop(temp) 
    }
  } else {
    # GLU format
    if (snp.list$stream == 0) {
      warning("Assuming snp.list$file.type is a GLU format. Setting snp.list$stream to 1")
      snp.list$stream <- 1
    }
  }

  # Check stream and nfiles
  if ((snp.list$stream) && (nfiles > 1)) {
    snp.list$stream <- 0
    stop("ERROR: snp.list$stream = 1 and length(snp.list$file) > 1")
  } 

  # For file type 1, stream must be 0
  if (snp.list$file.type == 1) {
    print("snp.list$stream is set to 0 for file.type = 1")
    snp.list$stream <- 0
  }
   
  # Check snpNames.list
  temp <- snp.list[["snpNames.list", exact=TRUE]]
  if (!is.null(temp)) {
    temp <- default.list(temp, c("file", "file.type", "delimiter", "header", "method"),
                         list("ERROR", 3, "\n", 0, 2), error=c(1, 0, 0, 0, 0))
    snp.list$snpNames.list <- temp
  }

  if (snp.list$file.type %in% c(9, 10)) snp.list$out.delimiter <- "\t"

  snp.list$alreadyChecked <- 1

  snp.list

} # END: check.snp.list

# Function to check pheno.list
check.pheno.list <- function(pheno.list) {

  if (is.null(pheno.list)) stop("ERROR: pheno.list must be specified")
  if (is.character(pheno.list)) pheno.list <- list(file=pheno.list, id.var="GWAS_ID", header=1)

  pheno.list <- default.list(pheno.list, 
   c("file", "id.var", "header", "remove.miss", "alreadyChecked", 
     "is.the.data", "in.miss", "orig.id.var", "unique.ids"),
        list("ERROR", "ERROR", 1, 0, 0, 0, c(NA, "NA", NaN, "NaN", "."), "origID5uh23gx25l6eq", 0), 
               error=c(1, 1, 0, 0, 0, 0, 0, 0, 0))

  if (pheno.list$alreadyChecked == 1) return(pheno.list)

  temp  <- getAllVars(pheno.list)
  pheno.list <- check.file.list(pheno.list, op=list(exist=1, vars=temp))

  if (!(pheno.list$file.type %in% c(1, 3, 4, 6, 8)))
    stop("ERROR: pheno.list$file.type must be 1, 3, 4, 6 or 8")

  # Check for sas.list
  if (pheno.list$file.type %in% c(4)) {
    if (is.null(pheno.list$sas.list)) stop("pheno.list$sas.list is not specified")
  }

  # Check for zip file with type 5
  if (pheno.list$file.type == 5) {
    if (is.null(pheno.list$zipFile)) stop("pheno.list$zipFile is not specified")
  }

  # Check id names list
  temp <- pheno.list[["keep.ids.list", exact=TRUE]]
  if (!is.null(temp)) {
    pheno.list$keep.ids <- getIdsFromFile(temp, 
                           id.vec=pheno.list[["keep.ids", exact=TRUE]])
    pheno.list$keep.ids.list <- NULL
  }

  # Check id names list
  temp <- pheno.list[["remove.ids.list", exact=TRUE]]
  if (!is.null(temp)) {
    pheno.list$remove.ids <- getIdsFromFile(temp, 
                           id.vec=pheno.list[["remove.ids", exact=TRUE]])
    pheno.list$remove.ids.list <- NULL
  }

  temp <- pheno.list[["keep.vars", exact=TRUE]]
  if (!is.null(temp)) {
    temp2 <- pheno.list[["factor.vars", exact=TRUE]]
    if (!is.null(temp2)) {
      if (!all(temp2 %in% temp)) stop("ERROR in check.pheno.list: all factor.vars not in keep.vars")
    }
  }

  pheno.list$alreadyChecked <- 1

  pheno.list

} # END: check.pheno.list

# Function to check locusMap.list
check.locusMap.list <- function(locusMap.list) {

  if (is.null(locusMap.list)) stop("ERROR: locusMap.list must be specified")

  locusMap.list <- default.list(locusMap.list, 
     c("file", "header", "snp.var","chrm.var", "loc.var", 
       "alreadyChecked"),
              list("ERROR", 0, "ERROR", "ERROR", "ERROR", 0), 
              error=c(1, 0, 1, 1, 1, 0))

  if (locusMap.list$alreadyChecked == 1) return(locusMap.list)

  locusMap.list$dir <- checkForSep(locusMap.list[["dir", exact=TRUE]])
  temp <- paste(locusMap.list$dir, locusMap.list$file, sep="")
  if (check.files(temp)) stop("ERROR: in check.locusMap.list") 

  if (is.null(locusMap.list[["file.type", exact=TRUE]])) {
    locusMap.list$file.type <- getFileType(temp[1], default=3)
  }
  if (is.null(locusMap.list[["delimiter", exact=TRUE]])) {
    locusMap.list$delimiter <- getFileDelim(temp[1], type=locusMap.list$file.type, default="")
  }


  if (!(locusMap.list$file.type %in% c(1, 3, 4, 6, 8)))
    stop("ERROR: locusMap.list$file.type must be 1, 3, 4, 6 or 8")

  if (!locusMap.list$header) {
    a1 <- check.vec.num(locusMap.list$snp.var, "locusMap.list$snp.var", len=1)
    a2 <- check.vec.num(locusMap.list$chrm.var, "locusMap.list$chrm.var", len=1)
    a3 <- check.vec.num(locusMap.list$loc.var, "locusMap.list$loc.var", len=1)
    if (sum(a1+a2+a3)) stop() 
  } else {
    locusMap.list$header <- 1
  }

  # Check for sas.list
  if (locusMap.list$file.type %in% c(4)) {
    if (is.null(locusMap.list$sas.list)) stop("locusMap.list$sas.list is not specified")
  }

  locusMap.list$alreadyChecked <- 1

  locusMap.list

} # END: check.locusMap.list

# Function to check temp.list
check.temp.list <- function(temp.list) {

  if (is.null(temp.list)) temp.list <- list()

  temp.list <- default.list(temp.list, c("dir", "delete", "id", "alreadyChecked"),
                               list(" ", 1, 1, 0))

  if (temp.list$alreadyChecked == 1) return(temp.list)

  # Check the directory
  temp.list$dir <- checkForSep(temp.list$dir) 
  if (temp.list$dir != "") {
    warn <- getOption("warn")
    # Change the warn option so a warning turns to an error
    options(warn=2)
    temp <- try(dir(temp.list$dir), silent=TRUE)
    if (class(temp) == "try-error") {
      temp <- paste("ERROR: the directory ", temp.list$dir, " does not exist", sep="")
      options(warn=warn)
      stop(temp)
    }
    options(warn=warn)
  }

  temp.list$alreadyChecked <- 1

  temp.list

} # END: check.temp.list

# Function to check if a file exists
check.files <- function(files) {
  # files   Character vector of files
  
  temp <- !(file.exists(files))
  if (any(temp)) {
    ret <- 1
    files <- files[temp]
    for (file in files) print(paste("The file ", file, " does not exist", sep=""))
  } else {
    ret <- 0
  }
  
  ret

} # END: check.files

# Function to check a numeric vector
check.vec.num <- function(vec, name, len=NULL, maxValue=NULL, 
                 minValue=NULL) {

  if (is.null(vec)) return(0)
  ret <- 0
  if (!is.numeric(vec)) {
    print(paste("ERROR: ", name, " must be of type numeric", sep=""))
    ret <- 1
  }

  if (!is.null(len)) {
    if (length(vec) != len) {
      print(paste("ERROR: ", name, " must be of length ", len, sep=""))
      ret <- 1
    }
  }

  if (!is.null(maxValue)) {
    if (max(vec) > maxValue) {
      print(paste("ERROR: each element of ", name, " must <= ", maxValue, sep=""))
      ret <- 1
    }
  }

  if (!is.null(minValue)) {
    if (min(vec) < minValue) {
      print(paste("ERROR: each element of ", name, " must >= ", minValue, sep=""))
      ret <- 1
    }
  }

  ret

} # END: check.vec.num

# Function to check a character vector
check.vec.char <- function(vec, name, len=NULL, checkList=NULL) {

  if (is.null(vec)) return(0)
  ret <- 0
  if (!is.character(vec)) {
    print(paste("ERROR: ", name, " must be of type character", sep=""))
    ret <- 1
  }

  if (!is.null(len)) {
    if (length(vec) != len) {
      print(paste("ERROR: ", name, " must be of length ", len, sep=""))
      ret <- 1
    }
  }

  if (!is.null(checkList)) {
    temp <- is.na(match(vec, checkList))
    if (sum(temp)) {
      ret  <- 1
      temp <- vec[temp]
      for (x in temp) {
        print(paste("ERROR: ", name, " contains the invalid element ", x, sep=""))
      }
    }
  }

  ret

} # END: check.vec.char

# Function to check for character vectors is a list
check.list.vec.char <- function(inList, names=NULL) {

  # inList    List to check
  # names     Character vector of names to check. If NULL, then
  #           all names are checked
  #           The default is NULL
  
  if (is.null(names)) names <- names(inList)
  sum <- 0
  for (name in names) {
    if (!is.null(inList[[name]])) {
      sum <- sum + check.vec.char(inList[[name]], name)
    }
  }
  if (sum) stop()
  0

} # END check.list.vec.char

# Function to check a vector
check.vec <- function(vec, name, opList) {

  opList <- default.list(opList, c("stopOnError"), list(1))

  if (is.numeric(vec)) {
    ret <- check.vec.num(vec, name, len=opList$len, 
            maxValue=opList$maxValue, minValue=opList$minValue)
  } else {
    ret <- check.vec.char(vec, name, len=opList$len, 
            checkList=opList$checkList)
  }
  if (ret) {
    if (opList$stopOnError) stop()
  }
  ret

} # END: check.vec

# Function to get the (approximate) number of rows to read
getRead.n <- function(file, what="character") {

  fid  <- file(file, "r")
  temp <- scan(file=fid, what="character", sep="\n", nlines=1, quiet=TRUE)
  p1   <- seek(fid, where=0, origin="start")  
  temp <- scan(file=fid, what="character", sep="\n", nlines=2, quiet=TRUE)
  p2   <- seek(fid, where=0, origin="start")
  p3   <- seek(fid, where=0, origin="end")
  p3   <- seek(fid, where=0, origin="end")
  close(fid)
  nrow <- p3/(p2-p1)
 
  

} # END: getRead.n

# Function to create dummy variables. The returned data object will
# have the original factor variables removed from it (see keep.factor)
# and new dummy variables added. For example if data[, "x1"] is a factor 
# with levels 0, 1, and 2, then the new variables will be "x1_1" and
# "x1_2", provided 0 is the baseline category.
# If you do not want to keep the baseline variable too, then set
#  baseline to a value that is not a level of the factor.
createDummy <- function(data, vars=NULL, baseline=NULL, 
                         keep.factor=NULL) {
  # data           Data frame, matrix, or vector
  # vars           Character vector of variable names or numeric vector
  #                of column numbers for  the factor variables that will be
  #                turned into dummy variables. If data is a vector, then
  #                vars does not have to be specified. 
  #                If NULL, then all the 
  #                factors will turn into dummy variables.
  #                The default is NULL.
  # baseline       Vector of baseline categories. The length of
  #                this vector must be equal to the length of vars or 
  #                the number of factors in data.
  #                The default is NULL.
  # keep.factor    Logical vector to keep the original factors in the 
  #                returned object. If NULL, then all factors are
  #                removed. The length of this vector must equal the
  #                length of vars or the number of factors. 
  #                The default is NULL.

  # Check for vector and column names
  if (is.null(dim(data))) {
    dim(data) <- c(length(data), 1)
    vars      <- "VAR1"
  }
  cnames <- colnames(data)
  nc     <- ncol(data)
  if (is.null(cnames)) {
    cnames <- paste("VAR", 1:nc, sep="")
    colnames(data) <- cnames
  }

  # Check for numeric vector vars
  if (is.numeric(vars)) vars <- cnames[vars]

  # Get the variables to factor
  if (is.null(vars)) {
    for (i in nc) {
      if (is.factor(data[, i])) vars <- c(vars, cnames[i])
    }
  } else {
    # Check the variable names
    temp <- check.vec.char(vars, "vars", len=NULL, checkList=colnames(data))
    if (temp) stop("ERROR in vars")
  }
  nvar <- length(vars)
  if (!nvar) return(list(data=data))

  if (!is.null(baseline)) {
    if (length(baseline) != nvar) stop("ERROR: baseline is incorrect")
    baseFlag <- 1
  } else {
    baseFlag <- 0
  }
  if (!is.null(keep.factor)) {
    if (length(keep.factor) != nvar) stop("ERROR: keep.factor is incorrect")
  } else {
    keep.factor <- rep(FALSE, times=nvar)
  }

  # Initialize a list
  newVars   <- list()
  ret       <- data
  nr        <- nrow(data)
  addedCols <- NULL

  for (i in 1:nvar) {
    var <- vars[i]

    # Get the levels
    temp <- data[, var]
    if (!is.factor(temp)) temp <- factor(temp)
    levels  <- levels(temp)
    if (length(levels) == 1) {
      keep.factor[i] <- TRUE
      next
    }

    # Get the baseline category
    if (!baseFlag) {
      base <- levels[1]
    } else {
      base <- baseline[i]
    }
     
    # Remove the baseline
    temp    <- !(levels == base) 
    levels  <- levels[temp]
    nlevels <- length(levels)
 
    # Initialize a dummy matrix
    temp           <- matrix(data=NA, nrow=nr, ncol=nlevels)
    addedVars      <- paste(var, "_", levels, sep="")
    colnames(temp) <- addedVars

    # Add the new variables names to the return list
    newVars[[var]] <- addedVars

    # Get the binary values
    for (j in 1:nlevels) {
      temp[, j] <- as.numeric(data[, var] == levels[j])
    }

    # Add the new variables
    addedCols <- cbind(addedCols, temp)  
  }

  # Add the new variables
  ret <- cbind(ret, addedCols)

  # Free memory
  rm(data, addedCols)
  temp <- gc(verbose=FALSE)

  # Remove original factors
  if (sum(keep.factor) != nvar) {
    temp <- vars[!keep.factor]
    ret  <- removeOrKeepCols(ret, temp, which=-1)
  }

  list(data=ret, newVars=newVars)

} # END: createDummy

# Function to match elements of 2 vectors
getOrder <- function(baseVec, newVec, removeMiss=0, errorIfMiss=1) {

  # baseVec        Baseline vector
  # newVec
  # removeMiss
  # errorIfMiss

  ret <- match(baseVec, newVec)
  if (errorIfMiss) {
    if (any(is.na(ret))) stop("ERROR: in matching vectors")
  }
  if (removeMiss) ret <- ret[!is.na(ret)]

  ret

} # END: getOrder

# Function to read in files and save as an R object file
dat2rda <- function(infile, outfile) {
                   
  dat <- readLines(infile)
  save(dat, file=outfile)

  0
} # END: dat2rda

# Function to read in file (numeric matrix) and save as an R object file
matrix2rda <- function(infile, outfile, delimiter="\t", stopOnError=1) {
                   
  tlist <- list(returnMatrix=1, include.row1=1, what=double(0),
                delimiter=delimiter, start.row=1, stop.row=-1)

  options(warn=2)
  dat <- try(scanFile(infile, tlist), silent=TRUE)
  options(warn=1)
  if (class(dat) == "try-error") {
    temp <- paste("ERROR with file ", infile, sep="")
    if (stopOnError) stop(temp)
    return(0)
  }
  save(dat, file=outfile)

  0
} # END: matrix2rda

# Function to get all ids requested
getIdsFromFile <- function(file.list, id.vec=NULL) {

  if (!is.null(file.list)) {
    file.list <- check.file.list(file.list)
    file.list <- default.list(file.list, c("id.var"), list(1))
    if (file.list$id.var == -1) {
      fid  <- getFID(file.list$file, file.list)
      temp <- scan(file=file.list$file, what="character", 
                   sep=file.list$delimiter)
      close(fid)
      id.vec <- c(id.vec, temp) 
    } else {
      var    <- file.list$id.var
      temp   <- loadData.table(file.list)
      id.vec <- c(id.vec, makeVector(temp[, var]))
    }
  }
  id.vec <- unique(id.vec)
  id.vec <- removeWhiteSpace(id.vec)
  if (!length(id.vec)) stop(paste("No ids in file", file.list$file))

  id.vec

} # END: getIdsFromFile

# Function to return the names of a vector or array
getNames <- function(obj) {

  if (is.null(dim(obj))) {
    ret <- names(obj)
  } else {
    ret <- dimnames(obj)
    ret <- ret[[1]]
  }
  ret

} # END: getNames

# Function to return or create variable names
getVarNames <- function(obj, prefix="VAR") {

  if (is.null(obj)) return(NULL)
  if (is.data.frame(obj)) return(colnames(obj))
  if (is.matrix(obj)) {
    ret <- colnames(obj)
    if (is.null(ret)) ret <- paste(prefix, 1:ncol(obj), sep="")
    return(ret)
  }
  ret <- names(obj)
  if (is.null(ret)) ret <- paste(prefix, 1:length(obj), sep="")
  ret

} # END: getVarNames

# Function to check if a directory ends with a slash.
# If not, then it will be added.
checkForSep <- function(dir) {

  if (is.null(dir)) return("")

  if (dir %in% c("", " ", "  ")) {
    dir <- paste(".", .Platform$file.sep, sep="")
    return("")
  }

  n <- nchar(dir)
  if (substring(dir, n, n) %in% c("\\", "/")) {
    return(dir)
  } else {
    dir <- paste(dir, .Platform$file.sep, sep="")
    return(dir)
  }

} # END: checkForSep

# Function to remove columns or variables from a matrix or data frame
removeOrKeepCols <- function(x, vars, which=1, stopOnError=1) {

  # x      Matrix or data frame 
  # vars   Vector of variable names or column numbers.
  # which  1 or -1 to keep or remove cols
  #        The default is 1 
  # stopOnError  0 or 1 to call stop() if columns are not found

  n      <- dim(x)
  dfFlag <- is.data.frame(x)
  if (is.null(n)) stop("ERROR: x should be 2 dimensional")
  if (which != 1) which <- -1

  if (!is.numeric(vars)) {
    cols <- match(vars, colnames(x))
    temp <- is.na(cols)
    if (any(temp)) {
      if (stopOnError) {
        print(vars[temp])
        print("The above columns were not found in the data")
        stop("ERROR in removeOrKeepCols")
      }
      cols <- cols[!temp]
      if (!length(cols)) {
        if (which == 1) {
          return(NULL)
        } else {
          return(x)
        }
      }
    }
  } else {
    cols <- vars
  }
  cols <- which*cols

  # Keep or remove columns
  ret <- x[, cols]
  
  # Check for NULL dimension
  if (is.null(dim(ret))) {
    dim(ret) <- c(n[1], length(ret)/n[1])
    if (dfFlag) ret <- data.frame(ret)
    cnames <- colnames(x)
    if (!is.null(cnames)) colnames(ret) <- cnames[cols]
    rownames(ret) <- rownames(x)
  }

  ret

} # END: removeOrKeepCols

# Function to remove columns or variables from a matrix or data frame
removeOrKeepRows <- function(x, rows, which=1) {

  # x      Matrix or data frame 
  # rows   Vector of row numbers, character vector of row names,
  #        or logical vector.
  # which  1 or -1 to keep or remove rows
  #        The default is 1 

  n      <- dim(x)
  dfFlag <- is.data.frame(x)
  if (is.null(n)) stop("ERROR: x should be 2 dimensional")
  if (which != 1) which <- -1
  rnames <- rownames(x)

  if (is.logical(rows)) {
    if (length(rows) != n[1]) stop("ERROR with logical vector rows")
    if (which != 1) rows <- !rows
  } else if (is.character(rows)) {
    rows <- match(rows, rnames)
    if (any(is.na(rows))) stop("ERROR in removeOrKeepRows")
    rows <- which*rows
  } else {
    rows <- which*rows
  }
  cnames <- colnames(x)

  # Keep or remove rows
  x <- x[rows, ]
  
  # Check for NULL dimension
  if (is.null(dim(x))) {
    dim(x) <- c(length(x)/n[2], n[2])
    if (dfFlag) x <- data.frame(x)
    if (!is.null(cnames)) colnames(x) <- cnames
    rownames(x) <- rnames[rows]
  }

  x

} # END: removeOrKeepRows

# Function to write snp rows to an open file
writeSnpLines <- function(snpNames, fid, snpData, nsub, orderFlag=0,
                    order=NULL, delimiter="|", sep="|") {

  # snpNames       Character vector of snps (no header)
  # fid            File connection
  # snpData        Character vector of the snp data
  # nsub           Number of subjects
  # orderFlag      0 or 1 if the subjects are already ordered
  #                0 = not ordered
  # order          Vector of integers for orderFlag = 1
  # delimiter      Input delimiter
  # sep            Output delimiter

  # Define a local function to search for character strings
  f1 <- function(str) {

    grep(str, snpData, value=FALSE)

  } # END: f1

  # Get the rows by searching for the snp names
  temp <- unlist(lapply(snpNames, f1))

  if (length(temp)) {
    snpData <- snpData[temp]
    nsnps   <- length(snpData)
  } else {
    return(0)    
  }

  nsubP1 <- nsub + 1
  ret    <- 0
  for (i in 1:nsnps) {
    temp  <- getVecFromStr(snpData[i], delimiter=delimiter)
    snp   <- temp[1]
    temp  <- temp[-1]
    if (!orderFlag) temp  <- temp[order]
    if (snp %in% snpNames) {
      write(c(snp, temp), file=fid, ncolumns=nsubP1, sep=sep)
      ret <- ret + 1
    }
  }

  ret

} # END: writeSnpLines

# Function to create factors in a data frame
factorVars <- function(data, vars) {

  # data    Data frame
  # vars    Vector of variables names or column numbers

  if (!is.data.frame(data)) stop("ERROR: data must be a data frame")
  if (is.numeric(vars)) vars <- colnames(data)[vars]

  for (var in vars) {
    data[, var] <- factor(data[, var])
  }

  data

} # END: factorVars

# Function to unfactor all columns ij a matrix or data frame
unfactor.all <- function(data, fun=NULL) {

  nc <- ncol(data)
  if (is.null(nc)) {
    nc <- 1
    dim(data) <- c(length(data), 1)
  }
  for (i in 1:nc) data[, i] <- unfactor(data[, i], fun=fun)
  data

} # END: unfactor.all

# Function to un-factor a factor
unfactor <- function(fac, fun=NULL) {

  # fac   Factor
  # fun   Function like as.character or as.numeric, etc

  if (is.factor(fac)) {
    ret <- levels(fac)[fac]
  } else {
    ret <- fac
  }

  if (!is.null(fun)) ret <- fun(ret)

  ret

} # END: unfactor

# Function to multiply each row or column of a matrix by a vector
matrixMultVec <- function(mat, vec, by=2) {

  # by    1 or 2  1 = rows, 2 = columns
  
  d <- dim(mat)
  if (by == 1) {
    vec <- rep(vec, each=d[1])
  } else {
    vec <- rep(vec, times=d[2])
  }

  dim(vec) <- d
  ret <- mat*vec
  ret

} # END: matrixMultVec

# Function to divide each row or column of a matrix by a vector
matrixDivideVec <- function(mat, vec, by=2) {

  # by    1 or 2  1 = rows, 2 = columns

  d <- dim(mat)
  if (by == 1) {
    vec <- rep(vec, each=d[1])
  } else {
    vec <- rep(vec, times=d[2])
  }

  dim(vec) <- d

  ret <- mat/vec
  ret

} # END: matrixDivideVec

# Function to get the column number of a 1 for each row in a matrix
#  of dummy variables. The matrix must contain only one 1 in each row.
# If the matrix has a 1 in the (i,j)th element, then in the returned
# vector ret, ret[i] = j.
getColNumber <- function(mat) {

  nc  <- ncol(mat)
  nr  <- nrow(mat)
  ret <- rep(NA, times=nr)
  for (i in 1:nc) {
    temp      <- mat[, i] == 1
    ret[temp] <- i
  }
  ret

} # END: getColNumber

# Function to return a block diagonal matrix
blockDiag <- function(mat.list) {
  
  # mat.list   List of matrices to form the block diagonal matrix

  # Get the dimensions of each matrix
  d <- matrix(unlist(lapply(mat.list, dim)), byrow=TRUE, ncol=2)

  # Initialize 
  ret  <- matrix(data=0, nrow=sum(d[,1]), ncol=sum(d[,2]))
  row0 <- 1
  row1 <- 0
  col0 <- 1
  col1 <- 0

  # Set each block
  for (i in 1:length(mat.list)) {
    row1 <- row1 + d[i,1]
    col1 <- col1 + d[i,2]
    ret[row0:row1, col0:col1] <- mat.list[[i]]
    row0 <- row1 + 1
    col0 <- col1 + 1
  }

  ret

} # END: blockDiag

# Function to add intercept column to a matrix or vector
addIntercept <- function(x, nrow=NULL) {

  if (is.null(x)) {
    if (is.null(nrow)) stop("Cannot add intercept")
    ret <- matrix(1, nrow=nrow, ncol=1)
  } else {
    if (is.null(nrow)) nrow <- nrow(x)
    if (is.null(nrow)) nrow <- length(x)
    nc <- ncol(x)
    if (is.null(nc)) {
      nc     <- 1
      dim(x) <- c(nrow, nc)
    } 
    ret <- cbind(rep(1, times=nrow), x)  
  }

  ret

} # END: addIntercept

# Function to use integers 1, ..., n for a categorical vector
changeStrata <- function(vec, start.n=1) {

  if (is.factor(vec)) vec <- unfactor(vec)
  uvec <- sort(unique(vec))
  ret  <- rep.int(0, times=length(vec))
  for (u in uvec) {
    ret[vec == u] <- start.n
    start.n <- start.n + 1
  }

  ret

} # END: changeStrata

# Function to remove rows that contain at least 1 missing value from
#  a data frame or matrix or vector
removeMiss <- function(x, miss=NA) {
  
  # x
  # miss    Vector of missing values

  d <- dim(x)

  # For a vector
  if (is.null(d)) {
    flag   <- 1
    cnames <- names(x)
    cFlag  <- !is.null(cnames)
    d      <- c(length(x), 1)
    dim(x) <- d
  } else {
    flag  <- 0
    cFlag <- 0
  }

  # Be careful for a data frame
  if (is.data.frame(x)) {
    temp <- matrix(data=TRUE, nrow=d[1], ncol=d[2])
    for (i in 1:d[2]) temp[, i] <- !(x[, i] %in% miss)
  } else {
    temp <- !(x %in% miss)
  }
  if (cFlag) cnames <- cnames[temp]

  dim(temp) <- d
  temp <- rowSums(temp)
  temp <- temp == d[2]
  x    <- removeOrKeepRows(x, temp, which=1)
  
  if (flag) {
    dim(x) <- NULL
    if (cFlag) names(x) <- cnames
  }
  
  x

} # END: removeMiss

# Function to remove missing values of a data frame from
#  certain variables
removeMiss.vars <- function(x, vars=NULL, miss=NA) {

  if (is.null(vars)) {
    x <- removeMiss(x, miss=miss)
    return(x)
  }

  temp <- rep(TRUE, times=nrow(x)) 
  for (var in vars) temp <- temp & !(x[, var] %in% miss)
  
  x <- removeOrKeepRows(x, temp, which=1)
  x

} # END: removeMiss.vars

# Function to write a (named) vector to a file. For type 3 and close = 0,
#  the file conection is returned.
writeVecToFile <- function(vec, file, colnames=NULL, type=3, close=1,
                           sep=" ") {

  # vec
  # file
  # colnames
  # type        1 or 3  
  # close       0 or 1 (for type 3 only)

  if (type == 1) {
    if (!is.null(colnames)) names(vec) <- colnames
    save(vec, file=file)
    return(0)
  } else {
    fid <- file(file, "w")
    n   <- length(vec)
    if (!is.null(colnames)) write(colnames, file=fid, ncolumns=n, sep=sep)
    write(vec, file=fid, ncolumns=n, sep=sep)
    if (close) {
      close(fid)
      fid <- 0
    }
    return(fid)
  }
  0

} # END: writeVecToFile 

# Function for writing vectors to files
writeVec <- function(vec, fileOrFID, colnames=NULL, isFID=0, sep="\t", 
                     close=1, type=3) {

  if (type == 1) {
    if (!is.null(colnames)) names(vec) <- colnames
    save(vec, file=fileOrFID)
    return(0)
  } 

  if (!isFID) fileOrFID <- file(fileOrFID, "w")
  if (!is.null(colnames)) {
    temp <- paste(colnames, collapse=sep, sep="")
    write(temp, file=fileOrFID, ncolumns=1)
  }
  temp <- paste(vec, collapse=sep, sep="")
  write(temp, file=fileOrFID, ncolumns=1)
  if (close) {
    close(fileOrFID)
    fileOrFID <- 0
  }

  fileOrFID

} # END: writeVec

# Function to return the vector of genotypes for the different modes
#  of inheritance
getInheritanceVec <- function(which, recode=1) {

  # which   NULL, 0, 1, 2, 3
  #         0 = trend
  #         1 = dominant
  #         2 = recessive
  #         3 = factor
  
  if ((!is.null(recode)) && (!recode)) return(NULL)

  if ((is.null(which)) || (!which)) {
    return(c(0, 1, 2))
  } else if (which == 1) {
    return(c(0, 1, 1))
  } else if (which == 2) {
    return(c(0, 0, 1))
  } else {
    return(c(0, 1, 2))
  }

} # END: getInheritanceVec

# Function to add interaction variables to a matrix or data frame
# Returns a list with the names data and newVars. newVars is a 
# character vector containing the variables added.
addInterVars <- function(data, vec, inter.vars, prefix="SNP") {

  # data           Data frame or matrix
  # vec            Numeric vector or factor for interactions
  #                If a matrix with more than 1 column, then
  #                it is assumed the columns are dummy variables
  # inter.vars     Matrix or data frame of variables that 
  #                will interact with var. These variables cannot
  #                be factors.
  # prefix         Variable prefix for interaction variables added
  #                The default is "SNP".

  cnames2 <- colnames(inter.vars) 
  if (is.null(cnames2)) cnames2 <- paste("V", 1:ncol(inter.vars), sep="")
 
  facFlag <- is.factor(vec)
  ncVec   <- ncol(vec)
  if (is.null(ncVec)) ncVec <- 0
  mFlag   <- ncVec > 1

  if (facFlag || mFlag) {
    newVars <- NULL
    if (facFlag) {
      vec <- data.frame(vec)
      colnames(vec) <- prefix
      vec <- createDummy(vec)$data
    }
    cnames <- colnames(vec)
    #temp2  <- NULL

    for (i in 1:ncol(vec)) {
      temp   <- matrixMultVec(inter.vars, vec[, i], by=2)
      tnames <- paste(cnames[i], "_", cnames2, sep="")
      newVars <- c(newVars, tnames) 
      colnames(temp) <- tnames
      if (i == 1) {
        temp2 <- temp
      } else {
        temp2 <- cbind(temp2, temp)
      }
    }
    data <- cbind(data, temp2)
  } else {
    temp <- matrixMultVec(inter.vars, vec, by=2)
    newVars <- paste(prefix, "_", cnames2, sep="")
    colnames(temp) <- newVars
    data <- cbind(data, temp)
  }
 
  list(data=data, newVars=newVars)

} # END: addInterVars

# Function to return variable names and positions from a vector
matchNames <- function(vec, name, exact=0) {

  # vec       Character vector to search from.
  # name      Character vector of names to search for.
  # exact     0 or 1 to perform exact matching
  #           The default is 0.

  names <- NULL
  pos   <- NULL

  if (exact) {
    temp <- match(name, vec)
    temp <- temp[!is.na(temp)]
    if (length(temp)) {
      names <- vec[temp]
      pos   <- temp
    } 
  } else {
    for (n in name) {
      temp <- grep(n, vec)
      if (length(temp)) {
        pos   <- c(pos, temp)
        names <- c(names, vec[temp])
      }
    }
    temp  <- !duplicated(pos)
    pos   <- pos[temp]
    names <- names[temp]
  }

  list(names=names, pos=pos)  

} # END: matchNames

# Function to get an option from the command arguments
getCommandArg <- function(optName, fun=NULL) {
  
  # optName      Unique option name in the command arguements
  # fun          Function to apply

  # Get the vector of command arguments
  comArgs <- commandArgs()

  # Get the one that we need
  temp <- grep(optName, comArgs)

  if (length(temp) == 1) { 
    ret <- sub(optName, "", comArgs[temp])
    if (!is.null(fun)) ret <- fun(ret)
  } else {
    temp <- paste("ERROR with ", optName, " in commandArgs()", sep="")
    stop(temp)
  }

  ret

} # END: getCommandArg

# Function to define list in swarm generator files
genfile.list <- function(inList, listName, fid) {

  # If the field is inList is a function or family,
  # put it in quotes and set the comment to "FUNCTION".

  names <- names(inList)
  llen  <- length(inList)
  if (is.null(names)) {
    flag  <- 1   
    names <- 1:llen
    str1  <- '[['
    str2  <- ']] <- '
  } else {
    flag <- 0
    str1 <- '[["'
    str2 <- '"]] <- '
  }
  temp  <- paste("\n ", listName, " <- list() \n", sep="")
  cat(temp, file=fid)

  if (!llen) return(NULL)
  for (name in names) {
    temp <- inList[[name]]
    cmm  <- comment(temp)
    if (is.null(cmm)) cmm <- "NULL"

    # For a list
    if (is.list(temp)) {
      if (flag) {
        name2 <- paste("list", name, sep="")
      } else {
        name2 <- name
      }
      genfile.list(temp, name2, fid)
      temp <- paste(listName, str1, name, str2, name2, ' \n', sep='')
      cat(temp, file=fid)
      next
    }
    if (length(temp) == 1) {
      if (is.character(temp)) { 
        if (cmm == "FUNCTION") {
          temp <- paste(listName, str1, name, str2, temp, ' \n', sep='')
        } else {
          temp <- paste(listName, str1, name, str2, '"', temp, '" \n', sep='')
        }
      } else {
        temp <- paste(listName, str1, name, str2, temp, ' \n', sep='')
      }
      cat(temp, file=fid) 
    } else {
      if (flag) {
        field <- paste(listName, "$", "'", name, "'", sep="")
      } else {
        field <- paste(listName, "$", name, sep="")
      }
      genfile.vec(temp, field, fid)
    }  
  }

  NULL

} # END: genfile.list

# Function to define a vector in swarm generator files
genfile.vec <- function(vec, name, fid) {

  n       <- length(vec)
  vecType <- "character"
  cflag   <- 1
  if (is.numeric(vec)) {
    vecType <- "numeric"
    cflag   <- 0
  }
  temp <- paste(name, " <- ", vecType, "(", n, ") \n", sep="")
  cat(temp, file=fid)
  k <- 1
  for (temp in vec) { 
    if (cflag) {
      temp <- paste(name, '[', k, '] <- "', temp, '" \n', sep='')
    } else {
      temp <- paste(name, '[', k, '] <- ', temp, ' \n', sep='')
    }
    cat(temp, file=fid)
    k <- k + 1
  }
  
  NULL

} # END: genfile.vec

# Function to sort a matrix or data frame by a column
sort2D <- function(data, col, dec=FALSE, fun=NULL) {

  # data   Data frame or matrix
  # col    Column number or variable name to sort on
  # dec    0 or 1 for decreasing
  #        The default is 0
  # fun    Function to apply to col before sorting
  #        The default is NULL

  vec  <- makeVector(unfactor(data[, col]))
  if (!is.null(fun)) vec <- fun(vec)
  temp <- is.na(vec)
  keep <- removeOrKeepRows(data, temp, which=1)
  data <- removeOrKeepRows(data, temp, which=-1)
  temp <- sort(vec, decreasing=dec, index.return=TRUE)$ix
  data <- rbind(data[temp, ], keep)
  return(data)

} # END: sort2D

# Function to update snp.list
update.snp.list <- function(snp.list, where=0) {

  # where    Integer specifying where in the program

  if (where == 1) {
    stream <- snp.list[["stream", exact=TRUE]]
    type   <- snp.list[["file.type", exact=TRUE]]
    if (stream) {
      if (type != 2) {
        n <- length(snp.list[["file", exact=TRUE]])
        snp.list$start.vec <- rep(1, times=n)
        snp.list$stop.vec  <- rep(-1, times=n)
      }
    }
  }

  temp <- snp.list[["snpNames.list", exact=TRUE]]
  if (!is.null(temp)) {
    snp.list$snpNames <- getIdsFromFile(temp, 
                           id.vec=snp.list[["snpNames", exact=TRUE]])
    snp.list$snpNames <- unique(snp.list$snpNames)
    snp.list$snpNames.list <- NULL
  }

  temp <- snp.list[["snpNames", exact=TRUE]]
  if (!is.null(temp)) {
    #n <- length(snp.list$file)
    #snp.list$start.vec <- rep(1, times=n)
    #snp.list$stop.vec  <- rep(-1, times=n)
  }

  snp.list

} # END: update.snp.list

# Function to return a temporary file name
getTempfile <- function(dir, prefix=NULL, ext=NULL) {

  pattern <- paste(prefix, collapse="", sep="")
  ret     <- tempfile(pattern=pattern, tmpdir="")
  ret     <- paste(dir, ret, ext, sep="") 

  ret

} # END: getTempfile

# Function to extract elements from a character vector
extractByStr <- function(dat, search, op=NULL) {

  # dat           Character vector to search in
  # search        Character vector for strings to search for in dat
  ################################################################
  # op            List with names:
  #  include.row1 0 or 1
  #               The default is 1
  #  substr.vec   Vector of max length 2 for start and stop options
  #               in the substring function
  #               The default is NULL.
  #  keep         0 or 1 to remove or keep matched strings
  #               The default is 1
  #  exact        0 or 1 for exact matching
  #               The default is 0
  #  delimiter    Delimiter in dat (used for exact matching)
  #               
  #  removeWhiteSpace 0 or 1 for exact matching
  #                   The default is 0

  # Define a local function to search for character strings
  f1 <- function(str) {

    grep(str, dat, value=FALSE)

  } # END: f1

  op <- default.list(op, c("include.row1", "keep", "exact", "removeWhiteSpace"), 
             list(1, 1, 0, 0))

  sep <- op[["delimiter", exact=TRUE]]
  if ((is.null(sep)) && (op$exact)) sep <- getFileDelim(dat) 
  keep <- op$keep

  sub.vec <- op[["substr.vec", exact=TRUE]]
  if (!is.null(sub.vec)) {
    subFlag <- 1
    a       <- sub.vec[1]
    b       <- sub.vec[2]
    save    <- dat
    if (is.na(b)) b <- Inf
    dat     <- substr(dat, a, b)
  } else {
    subFlag <- 0
  }

  # Get the rows by searching for the snp names
  rows <- unlist(lapply(search, f1))
  nr   <- length(rows)
  if (!nr) return(NULL)
 
  # Get a logical vector of rows to keep or drop
  temp <- (1:length(dat)) %in% rows
  n    <- length(temp)

  # For exact matching
  if ((op$exact) && (n > 1)) {
    id <- character(n)
    for (i in 2:n) {
      if (temp[i]) id[i] <- getVecFromStr(dat[i], delimiter=sep)[1] 
    }
    if (op$removeWhiteSpace) {
      search <- removeWhiteSpace(search)
      id     <- removeWhiteSpace(id)
    }
    temp <- temp & (id %in% search)
    rm(id)
    gc()
  }

  if (keep == 0) temp <- !temp 

  # Add the first row if needed
  if (op$include.row1) {
    temp[1] <- TRUE
  } else {
    temp[1] <- FALSE
  }

  # Get the subset
  if (subFlag) {
    return(save[temp])
  } else {
    return(dat[temp])
  }

} # END: extractByStr

# Function to close and delete a file
closeFile <- function(fid, file=NULL, delete=0) {

  close(fid)
  ret <- 0
  if ((delete) && (!is.null(file))) {
    if (file.exists(file)) ret <- file.remove(file)
  }
  ret

} # END: closeFile

# Function to create a vector from a matrix
makeVector <- function(x) {

  d <- dim(x)
  if (is.null(d)) return(x)

  nn <- NULL
  if (d[1] == 1) {
    nn <- colnames(x)
  } else if (d[2] == 1) {
    nn <- rownames(x)
  }
  dim(x) <- NULL
  if (!is.null(nn)) names(x) <- nn
  if ((!is.vector(x)) && (is.list(x))) x <- unlist(x)
  x
 
} # END: makeVector

# Function to initialize a data frame
initDataFrame <- function(nrow, columns, rownames=NULL, colnames=NULL,
                  initChar="", initNum=NA) {
  
  cc   <- character(nrow)
  nn   <- double(nrow)
  nn[] <- initNum
  cc[] <- initChar
  cc   <- data.frame(cc)
  nn   <- data.frame(nn)
  nc   <- length(columns)

  columns <- toupper(columns)
  if (columns[1] == "C") {
    data <- cc
  } else {
    data <- nn
  }
  
  if (nc > 1) { 
    for (x in columns[-1]) {
      if (x == "C") {
        data <- cbind(data, cc)
      } else {
        data <- cbind(data, nn)
      }
    }
  }  
  data <- data.frame(data)
  if (!is.null(rownames)) rownames(data) <- rownames
  if (!is.null(colnames)) colnames(data) <- colnames

  for (i in 1:nc) {
    if (is.factor(data[, i])) data[, i] <- unfactor(data[, i])
  }
  data

} # END: initDataFrame

# Function to subset data by a single variable
subsetData.var <- function(data, var, operator, value, which=1, 
                           returnRows=0, na.value=FALSE) {

  # data
  # var
  # operator
  # value
  # which        1 or -1 to keep or drop
  #              The default is 1
  # returnRows   0 or 1
  #              The default is 0
  # na.value     TRUE or FALSE on what to do with NAs
  #              The default is FALSE

  lenv <- length(value)
  if (any(is.na(value))) {
    if (lenv > 1) stop("ERROR: IN subsetData.var: with NA in value")
    naFlag  <- 1
    numFlag <- 0
  } else {
    numFlag <- is.numeric(value)
    naFlag  <- 0
  }

  # Be careful with a vector for value
  if (lenv > 1) {
    if (numFlag) {
      value <- paste(value, collapse=",", sep="")
      value <- paste("c(", value, ")", sep="")
    } else {
      # Change the operator
      if (operator == "%in%") {
        operator <- "=="
      } else { 
        stop("ERROR: value cannot be a character vector")
      }
    }
  }

  if (naFlag) {
    temp <- ""
    if (operator == "!=") temp <- "!"  
    callStr <- paste("(", temp, "is.na(data[, var]))", sep=" ") 
  }
  else if (numFlag) {
    callStr <- paste("(as.numeric(data[, var])", operator, value, ")", sep=" ") 
  } else {
    callStr <- ""
    for (i in 1:lenv) {
      temp <- paste('(data[, var] ', operator, ' "', value[i], '")', sep='')
      if (i < lenv) temp <- paste(temp, " | ", sep="")
      callStr <- paste(callStr, temp, sep="")
    }
  }

  rows <- eval(parse(text=callStr))
  rows[is.na(rows)] <- na.value
  if (returnRows) return(rows)

  data <- removeOrKeepRows(data, rows, which=which)
  data

} # END: subsetData.var

# Function to subset data. Returns the subsetted data, or if op$which = 0,
#   returns the vector of TRUE/FALSE to subset
subsetData.list <- function(data, slist, returnRows=0) {

  # data
  ###################################################################
  # slist       List of sublists with names:
  #  var
  #  operator
  #  value
  #  which      1 or -1 to include or NOT within each list
  #             The default is 1
  #  logic.op   Logical operator 
  #             Only used starting from the second sublist
  #             The default is "&"
  #  na.value   TRUE or FALSE
  #             The default is FALSE
  #  last.which 1 or -1 If -1, the rows defined by all the sublists will
  #             be negated. The value for last.which in the last sublist
  #             is the only one applied.
  #             The default is 1.
  ####################################################################
  # returnRows  Set to 1 to return the logical vector of rows instead of
  #             the data.
  #             The default is 0

  n      <- length(slist)
  wvec   <- rep(9999, times=n)
  cnames <- colnames(data)
  cflag  <- !is.null(cnames)

  for (i in 1:n) {

    tlist    <- slist[[i]]
    if (!is.list(tlist)) stop("ERROR in subsetData.list: Input list is incorrect")
    tlist    <- default.list(tlist, 
                  c("var", "operator", "value", "which", "logic.op", "na.value", "last.which"), 
                  list("ERROR", "ERROR", "ERROR", 1, "&", FALSE, 1), error=c(1,1,1,0,0,0,0))
    var      <- tlist[["var", exact=TRUE]]
    if ((cflag) && (is.character(var))) {
      if (!(var %in% cnames)) {
        temp <- paste("ERROR in subsetData.list: ", var, " not in data", sep="")
        print(temp)
        stop()
      }
    }
    operator <- tlist[["operator", exact=TRUE]]
    value    <- tlist[["value", exact=TRUE]]
    which    <- tlist[["which", exact=TRUE]]
    last     <- tlist[["last.which", exact=TRUE]]
    wvec[i]  <- which
    na.value <- tlist[["na.value", exact=TRUE]]
    temp     <- subsetData.var(data, var, operator, value, returnRows=1, na.value=na.value)
    if (which == -1) temp <- !temp
    if (i == 1) {
      rows <- temp
    } else {
      logic.op <- tlist[["logic.op", exact=TRUE]]
      callStr  <- paste("rows ", logic.op, " temp", sep="") 
      rows     <- eval(parse(text=callStr))      
    }
  }

  if (last == -1) rows <- !rows
  if (returnRows) return(rows)

  data <- removeOrKeepRows(data, rows, which=1)

  data

} # END: subsetData.list

# Function to rename a variable on a matrix or data frame
renameVar <- function(data, old.var, new.var) {

  cnames <- colnames(data)
  i      <- match(old.var, cnames)
  if (is.na(i)) return(data)
  cnames[i] <- new.var
  colnames(data) <- cnames
  data

} # renameVar

# Function to return a list of variable names for a SNP
getVarNames.snp <- function(prefix="SNP_", genetic.model=0) {

  # genetic.model = 3 is for factors, assuming 0 is the baseline
  if (genetic.model != 3) {
    return(prefix)
  } else {
    return(paste(prefix, 1:2, sep=""))
  }

} # END: getVarNames.snp

# Function to return a list of variable names for interactions
getVarNames.int <- function(V, prefix="SNP_", genetic.model=0, sep="_") {

  # V   Matrix of interactions

  vnames <- getVarNames(V, prefix="V")

  # genetic.model = 3 is for factors, assuming 0 is the baseline
  if (genetic.model != 3) return(paste(prefix, sep, vnames, sep="")) 
  
  ret <- c(paste(prefix, 1, sep, vnames, sep=""),
           paste(prefix, 2, sep, vnames, sep=""))

  ret

} # END: getVarNames.int

# Function to check for a constant variable in a matrix or data frame
# Returns a list of the data with the constant variables removed, and
#  a vector of variable names removed
checkForConstantVar <- function(data, msg=1, removeVars=1) {

  # data
  # msg          0, 1, 2  0 = no message, 1 = warning, 2 = error
  # removeVars   0 or 1 to remove the constant variables from the
  #              data.
  #              The default is 1

  temp   <- apply(data, 2, "var", na.rm=TRUE)
  remove <- NULL
  temp2  <- temp == 0
  temp2[is.na(temp2)] <- FALSE
  ret    <- (1:length(temp))[temp2]

  if (length(ret)) {
    temp <- colnames(data)
    if (!is.null(temp)) {
      temp <- temp[ret]
    } else {
      temp <- ret
    }
    remove <- temp
    if (removeVars) data <- removeOrKeepCols(data, remove, which=-1)    

    # Message
    temp <- paste(temp, collapse=",")
    if (msg == 1) {
      temp <- paste("WARNING: Variables ", temp, " are constant and have been removed from the data", sep="") 
      warning(temp)
    } else if (msg == 2) {
      temp <- paste("ERROR: Variables ", temp, " are constant", sep="") 
      stop(temp)
    }
  } 

  list(data=data, remove=remove) 

} # END: checkForConstantVar

# Function to return a logical matrix from a stratification variable
strataMatrix <- function(strata) {

  # strata    Stratification vector

  strata <- unfactor(strata)
  nr     <- length(strata)
  us     <- unique(strata)
  nc     <- length(us)
  ret    <- matrix(data=FALSE, nrow=nr, ncol=nc)
  colnames(ret) <- us
  
  for (i in 1:nc) ret[, i] <- (strata == us[i])
  
  ret

} # END: strataMat

rep.cols <- function(vec, times) {
  # Function to replicate columns
  len      <- length(vec)
  dim(vec) <- c(len, 1)
  mat      <- rep(vec, times)
  dim(mat) <- c(len, times)
  mat
} # END: rep.cols

rep.mat <- function(scalar, nrow=NULL, ncol=NULL)
{
  mat      <- rep(scalar, each=nrow, times=ncol)
  dim(mat) <- c(nrow, ncol)
  mat
}

rep.rows <- function(vec, times) {
  # Function to replicate rows
  len      <- length(vec)
  dim(vec) <- c(1, len)
  mat      <- rep(vec, each=times)
  dim(mat) <- c(times, len)
  mat
} # END: rep.rows

# Function to merge 2 data frames 
mergeData <- function(base.data, new.data, new.vars, 
                      base.id="id", new.id="id") {

  # base.data   Data frame
  # new.data    Data frame
  # new.vars
  # base.id
  # new.id

  new.vars <- unique(new.vars)
  base.id  <- unique(base.id)
  new.id   <- unique(new.id)
  n.base   <- length(base.id)
  n.new    <- length(new.id)
  flag     <- 0

  if (n.base != n.new) stop("Number of id variables do not match")
  #if (nrow(new.data) != length(unique(new.data[, new.id]))) {
  #  stop("Ids not unique in new.data")
  #}
  if (nrow(base.data) != length(unique(base.data[, base.id]))) {
    flag <- 1
  }
  if (n.base > 1) flag <- 1

  temp <- match(new.id, new.vars)
  if (!is.na(temp)) new.vars <- new.vars[-temp]
  temp <- new.vars %in% colnames(base.data)
  if (any(temp)) {
     new.vars2 <- c(new.vars[!temp], paste(new.vars[temp], "_2", sep=""))
  } else {
     new.vars2 <- new.vars
  }

  # Add variables
  for (var in new.vars2) base.data[, var] <- NA
  nvars <- length(new.vars2)

  if (!flag) {
    rows <- base.data[, base.id] %in% new.data[, new.id]
    temp <- match(new.data[, new.id], base.data[, base.id])
    temp <- temp[!is.na(temp)]
    for (i in 1:nvars) {
      base.data[rows, new.vars2[i]] <- new.data[temp, new.vars[i]]
    } 
  } else {
    n <- nrow(base.data)
    for (i in 1:nrow(new.data)) {
      temp <- rep.int(TRUE, times=n)
      for (j in 1:n.base) {
        temp <- temp & (base.data[, base.id[j]] == new.data[i, new.id[j]])
      }
      if (any(temp)) base.data[temp, new.vars2] <- new.data[i, new.vars]
    }
  }

  base.data

} # END: mergeData

# Function to call the operating system
callOS <- function(command, intern=FALSE) {

  # Determine the platform
  os      <- .Platform$OS.type
  winFlag <- (os == "windows")

  if (winFlag) {
    ret <- shell(command, intern=intern)
  } else {
    ret <- system(command, intern=intern)
  }
  ret

} # END: callOS

# Function to get columns from a character vector
getColsFromCharVec <- function(vec, cols, delimiter="\t", colNames=NULL,
                        fun=as.numeric, ret.type=3) {

  # vec
  # cols        Numeric or character vector of columns 
  # colNames    NULL or ordered vector of column names in vec
  #             if cols is character
  #             The default is NULL

  cnamesFlag <- !is.null(colNames)

  if (is.character(cols)) {
    cols <- match(cols, colNames)
  }
  n <- length(vec)

  if (ret.type == 2) {
    ret <- vec
    for (i in 1:n) {
      temp   <- getVecFromStr(vec[i], delimiter=delimiter)
      ret[i] <- paste(temp[cols], collapse=delimiter, sep="")
    }
  } else {
    if (n == 1) {
      temp <- getVecFromStr(vec, delimiter=delimiter)
      ret  <- fun(temp[cols])
      if (cnamesFlag) names(ret) <- colNames[cols] 
    } else {
      ret <- matrix(data=NA, nrow=n, ncol=length(cols))
      if (cnamesFlag) colnames(ret) <- colNames[cols] 
      for (i in 1:n) {
        temp     <- getVecFromStr(vec[i], delimiter=delimiter)
        ret[i, ] <- fun(temp[cols])
      }
    }
  }

  ret

} # END: getColsFromCharVec 

# Function to search and replace strings within the names of an oject
changeStr.names <- function(obj, search, replace="") {

  # obj
  # search
  # replace

  d <- dim(obj)
  if (is.null(d)) {
    n <- names(obj)
    if (!is.null(n)) names(obj) <- gsub(search, replace, n, fixed=TRUE)
  } else {
    n <- colnames(obj)
    if (!is.null(n)) colnames(obj) <- gsub(search, replace, n, fixed=TRUE)
    n <- rownames(obj)
    if (!is.null(n)) rownames(obj) <- gsub(search, replace, n, fixed=TRUE)
  }
  obj

} # END: changeStr.names

# Function to get the correct phenotype data when there are formulas
#  to be applied. The original phenotype data will be returned with
applyFormulas <- function(data, formulas, remove=NULL) {

  # data       Data frame
  # formulas   List of formulas with variables in the data frame

  flag   <- 0
  rflag  <- !is.null(remove)
  data2  <- data
  ids    <- 1:nrow(data2)
  rownames(data2) <- ids 
  for (i in 1:length(formulas)) {
    f <- formulas[[i]]
    if ("formula" %in% class(f)) {
      # Get the design matrix
      temp <- model.matrix(f, data=data2)
      ids  <- intersect(ids, rownames(temp))
      flag <- 1
      if (rflag) {
        temp <- removeMiss(data.frame(temp), miss=remove)
        ids  <- intersect(ids, rownames(temp))
      }
    }
  } 

  if (flag) data <- removeOrKeepRows(data, as.numeric(ids), which=1)
  data

} # END: applyFormulas

# Function to return the formulas from a list
getFormulas <- function(inlist) {

  # inlist    List

  ret   <- list()
  index <- 1
  for (i in 1:length(inlist)) {
    if ("formula" %in% class(inlist[[i]])) {
      ret[[index]] <- inlist[[i]]
      index        <- index + 1
    }
  }
  ret

} # END: getFormulas

# Function to return the variables from a particular object
getAllVars <- function(obj, names=c("response.var", "main.vars", "int.vars",
                                    "strata.var", "start.var", "stop.var",
                                    "partition.var", "group.var", "id.var",
                                    "keep.vars", "remove.vars", "factor.vars")) {

  if (is.null(obj)) return(NULL)
  clss <- class(obj)

  # Character vector
  if ((is.vector(obj)) && ("character" %in% clss)) return(unique(obj))

  # Formula
  if ("formula" %in% clss) return(unique(all.vars(obj)))

  # List
  if ("list" %in% clss) {
    ret <- NULL
    if (is.null(names)) names <- 1:length(obj)
    for (nn in names) {
      obj.n <- obj[[nn, exact=TRUE]]
      if (is.null(obj.n)) next

      clss <- class(obj.n)      
      if ((is.vector(obj)) && ("character" %in% clss)) {
        ret <- c(ret, obj.n)
      } else if ("formula" %in% clss) {
        ret <- c(ret, all.vars(obj.n))
      }
    }
    slist <- obj[["subsetData", exact=TRUE]]
    if (!is.null(slist)) ret <- c(ret, getSubsetDataVars(slist))

    return(unique(ret))
  }

  stop("ERROR in formulaVars: obj is of wrong type")

} # END: getAllVars

# Function to change the alleles in a data frame
changeAlleles <- function(data, alleleVars, rows=NULL, newVars=NULL) {

  # data
  # alleleVars
  # rows         Logical vector of length = nrow(data) or NULL
  #              The default is NULL, so that all rows will be changed
  # newVars      New variable names. If NULL, then allelVars will be used
  #              The order of newVars must match alleleVars.
  #              The default is NULL.

  if (is.null(rows)) rows <- rep(TRUE, times=nrow(data))
  if (is.null(newVars)) newVars <- alleleVars
  i <- 1
  for (var in alleleVars) {
    tempA <- (data[, var] == "A") & rows
    tempT <- (data[, var] == "T") & rows
    tempG <- (data[, var] == "G") & rows
    tempC <- (data[, var] == "C") & rows
    nv    <- newVars[i]
    data[tempA, nv] <- "T"
    data[tempT, nv] <- "A"
    data[tempG, nv] <- "C"
    data[tempC, nv] <- "G"
    i <- i + 1
  }
  data

} # END: changeAlleles

# Function to change the genotypes in a data frame
changeGenotypes <- function(data, genoVars, rows=NULL, newVars=NULL) {

  # data
  # genoVars     Must be coded as "AA", "GG", "CT" etc
  # rows         Logical vector of length = nrow(data) or NULL
  #              The default is NULL, so that all rows will be changed
  # newVars      New variable names. If NULL, then allelVars will be used
  #              The order of newVars must match alleleVars.
  #              The default is NULL.

  if (is.null(rows)) rows <- rep(TRUE, times=nrow(data))
  if (is.null(newVars)) newVars <- genoVars

  all  <- c("AA", "CC", "GG", "TT", 
           "AC", "CA", "AG", "GA", "AT", "TA",
           "CG", "GC", "CT", "TC",
           "GT", "TG")
  all2 <- c("TT", "GG", "CC", "AA", 
           "TG", "GT", "TC", "CT", "TA", "AT",
           "GC", "CG", "GA", "AG",
           "CA", "AC")

  index <- 1
  for (var in genoVars) {
    x     <- makeVector(data[, var])
    genos <- unique(x[rows])
    ngeno <- length(genos)
    temp  <- all %in% genos
    if (sum(temp) != ngeno) next
    a0    <- all[temp]
    a2    <- all2[temp]

    nv    <- newVars[index]
    for (i in 1:ngeno) {
      temp <- (x %in% a0[i]) & rows
      data[temp, nv] <- a2[i]
    }
    index <- index + 1
   
  }
  data

} # END: changeGenotypes

# Function to write a table
writeTable <- function(x, outfile, delimiter="\t") {

  write.table(x, file=outfile, sep=delimiter, row.names=FALSE, quote=FALSE)

} # END: writeTable

# Function to compute cross tab for ids data in the analysis
crossTab <- function(snp.list, pheno.list, varlist, op=NULL) {

  # snp.list
  # pheno.list
  # varlist       List of character vectors of length 2 
  ###########################################################################
  # op            List with names
  #  outfile
  #  by.var       By variable in pheno.list data
  #  exclude      Character string of values to exclude
  #  temp.list

  snpFlag <- !is.null(snp.list)
  if ((!snpFlag) && (is.null(pheno.list[["id.var", exact=TRUE]]))) {
    pheno.list$id.var <- 1
  }

  # Check the input lists
  pheno.list <- check.pheno.list(pheno.list)

  outfile <- op[["outfile", exact=TRUE]]
  outflag <- !is.null(outfile)
  by.var  <- op[["by.var", exact=TRUE]]
  byflag  <- !is.null(by.var)

  if (snpFlag) {
    snp.list  <- check.snp.list(snp.list)
    temp.list <- op[["temp.list", exact=TRUE]]
  
    pheno.list$remove.miss <- 0
    pheno.list$make.dummy  <- 0
    temp <- pheno.list[["keep.vars", exact=TRUE]]
    if (!is.null(temp)) {
      pheno.list$keep.vars <- c(pheno.list$id.var, temp, by.var)
    }
    snp.list$stop.vec <- 2

    # Get the data vector of snps
    tlist <- list(include.row1=0, include.snps=0, return.type=1,
                missing=1, snpNames=1, orderByPheno=1, return.pheno=1)

    temp  <- try(getData.1(snp.list, pheno.list, temp.list, op=tlist),
               silent=TRUE)
    if (class(temp) == "try-error") {
      print(temp)
      stop("ERROR loading data")
    }

    # Get the phenotype data
    phenoData.list <- temp$phenoData.list
    x              <- phenoData.list$data

    rm(temp, phenoData.list, tlist)
    temp <- gc()
  } else {
    x <- loadData(pheno.list$file, pheno.list)
    pheno.list$data        <- data
    pheno.list$is.the.data <- 1
    x <- getPhenoData(pheno.list)$data 
    rm(pheno.list)
    gc()
  }

  if (outflag) sink(outfile)
  if (byflag) {
    byLevels <- unique(unfactor(x[, by.var]))
  } else {
    byLevels <- 1
  }

  len <- length(varlist)
  n   <- nrow(x)
  exclude <- op[["exclude", exact=TRUE]]
  if (is.null(exclude)) exclude <- "NULL"
  suffix <- paste(", exclude=", exclude, ")", sep="") 

  for (level in byLevels) {
    if (byflag) {
      temp <- (x[, by.var] == level)
      print(level)
    } else {
      temp <- rep(TRUE, times=n)
    }
    
    for (i in 1:len) {
      vars <- varlist[[i]]
      str  <- parse.vec(length(vars), vec.prefix="x[temp, vars[", vec.suffix="]]",
                   prefix="table(", suffix=suffix, delimiter=", ")
      tab  <- eval(parse(text=str))
      print(vars)
      print(tab) 
    }
  }
  
  if (outflag) sink()

  0

} # END: crossTab

# Function to merge the phenotype and genotype data
mergePhenoGeno <- function(snp.list, pheno.list, temp.list=NULL, op=NULL) {

  # Returns a list with names data and geno.vars

  # snp.list
  # pheno.list
  # temp.list
  # op
  #   outfile
  #   which    Vector of 0-2.  0=expected genotype count, 1=original coding, 
  #            2=dummy vars with NAs in the dummy vars
  #            Then default is 0
  #   alleles  (For snp.list$file.type = 9 or 10)
  #            0 or 1 to add alleles (major/minor) for each SNP
  #            The default is 0

  op <- default.list(op, c("which", "alleles"), list(0, 0))

  # Check the input lists
  snp.list   <- check.snp.list(snp.list)
  pheno.list <- check.pheno.list(pheno.list)

  imputeFlag <- snp.list$file.type %in% c(9, 10)
  snp.list$recode <- 0
  snp.list$genetic.model <- 3
  flag0 <- 0 %in% op$which
  flag1 <- 1 %in% op$which
  flag2 <- 2 %in% op$which
  aflag <- 0
  if ((imputeFlag) && (op$alleles)) aflag <- 1

  # Get the data vector of snps
  tlist <- list(include.row1=0, include.snps=0, return.type=1,
                missing=1, snpNames=1, orderByPheno=1, return.pheno=1)

  temp  <- try(getData.1(snp.list, pheno.list, temp.list, op=tlist),
               silent=TRUE)
  if (class(temp) == "try-error") {
    print(temp)
    stop("ERROR loading data")
  }

  snpData   <- temp$data
  snpNames  <- temp$snpNames
  delimiter <- getDelimiter(snp.list, output=1)
  nsnps     <- length(snpData)
  if (aflag) alleles <- temp$alleles

  # Get the phenotype data
  phenoData.list <- temp$phenoData.list
  phenoData0     <- phenoData.list$data

  # Determine if cc.var was specified
  ccVar <- pheno.list[["cc.var", exact=TRUE]]
  if (!is.null(ccVar)) {
    subset <- phenoData0[, ccVar] == 0
    subset[is.na(subset)] <- FALSE
  } else {
    subset <- rep(TRUE, times=nrow(phenoData0))
  }

  hcodes <- snp.list$heter.codes
  nc     <- ncol(phenoData0)

  for (i in 1:nsnps) {
    snp.v <- snpNames[i]
    snp   <- getVecFromStr(snpData[i], delimiter=delimiter)

    if ((!imputeFlag) && (any(op$which %in% c(0, 2)))) {
      exsnp <- as.integer(recode.geno(snp, in.miss="NA", subset=subset, heter.codes=hcodes)$vec)
    } else if (imputeFlag) {
      temp <- as.numeric(unlist(strsplit(snp, "|", fixed=TRUE)))
      mat  <- matrix(temp, byrow=TRUE, ncol=3)
    }

    if (flag0) {
      # Expected genotype counts
      if (!imputeFlag) {
        temp <- exsnp
      } else {
        temp <- mat[, 2] + 2*mat[, 3]
      }
      phenoData0[, snp.v] <- temp
    }
    if (flag1) {
      # Original
      new <- paste(snp.v, "_GENO", sep="")
      phenoData0[, new] <- snp
    }

    if (flag2) {
      if (imputeFlag) {
        phenoData0[, paste(snp.v, "_0", sep="")] <- mat[, 1]
        phenoData0[, paste(snp.v, "_1", sep="")] <- mat[, 2]
        phenoData0[, paste(snp.v, "_2", sep="")] <- mat[, 3]
      } else {
        phenoData0[, paste(snp.v, "_0", sep="")] <- as.numeric(exsnp == 0)
        phenoData0[, paste(snp.v, "_1", sep="")] <- as.numeric(exsnp == 1)
        phenoData0[, paste(snp.v, "_2", sep="")] <- as.numeric(exsnp == 2)
      }
    }
    if (aflag) {
      new <- paste(snp.v, "_ALLELES", sep="")
      phenoData0[, new] <- alleles[i]
    }
  }

  out <- op[["outfile", exact=TRUE]]
  if (!is.null(out)) writeTable(phenoData0, out)

  newnc  <- ncol(phenoData0)
  gnames <- colnames(phenoData0)[(nc+1):newnc]

  list(data=phenoData0, geno.vars=gnames)

} # END: mergePhenoGeno

# Function to parse a vector. Returns the string to evaluate.
parse.vec <- function(vec.len, vec.prefix="", vec.suffix="",
                      prefix="", suffix="", delimiter=", ") {

  # vec.len
  # vec.prefix    Constant string (for now)
  # vec.suffix    Constant string (for now)

  str <- prefix 
  for (j in 1:vec.len) {
    str <- paste(str, vec.prefix, j, vec.suffix, sep="")
    if (j < vec.len) str <- paste(str, delimiter, sep="")
  }
  str <- paste(str, suffix, sep="")
  str

} # END: parse.vec

# Function to replace strings in levels of a categorical variable
replaceStr.var <- function(data, var, str=" ", newStr="") {

  data[, var] <- gsub(str, newStr, data[, var], fixed=TRUE)

  data

} # END: replaceStr.var

# Function to replace strings in a data frame
replaceStr.list <- function(data, varlist) {

  for (i in 1:length(varlist)) {
    tlist <- varlist[[i]]
    tlist <- default.list(tlist, c("var", "str", "newStr"),
                          list("ERROR", " ", ""), error=c(1, 0, 0))
    data[, tlist$var] <- gsub(tlist$str, tlist$newStr, data[, tlist$var], fixed=TRUE)
  }
  data

} # END: replaceStr.var

# Function to add columns to a data frame from another file or data frame
# A data frame will be returned
addColumn <- function(data, id.var, file.list, op=NULL) {

  # data        Data frame or matrix. 
  # id.var      ID variables on data (either 1 or 2 variables)
  #             No default
  # file.list   List of type file.list with the additional names
  #   id.var    Id variables to match against id.var
  #             The order of this vector must match id.var.
  #             No default
  #   vars      Vars to copy to data.
  #             No default
  #   names     New variable names of vars to copy to data
  #             The default is vars
  #   type      Vector of "C" or "N" for character/numeric variables
  #             The default is NULL
  #######################################################################
  # op          List with names:
  #  leading0   0 or 1 to check for leading zeros in the id variable
  #             The default is 0
  #  initValue  Initial value for new columns
  #             The default is NA.
  #  replace    0 or 1 to replace the new variables
  #             The default is 1

  op <- default.list(op, c("leading0", "initValue", "replace"), list(0, NA, 1))

  file.list <- default.list(file.list,
               c("file", "id.var", "vars"), 
               list("ERROR", "ERROR", "ERROR"), 
               error=c(1, 1, 1))
  temp <- list(vars=c(file.list$id.var, file.list$vars))
  file.list <- check.file.list(file.list, op=temp)

  check.vec(id.var, "id.var", list(checkList=colnames(data), maxValue=ncol(data), minValue=1))

  type  <- file.list[["type", exact=TRUE]]
  vars  <- file.list$vars
  id    <- file.list$id.var
  names <- file.list$names
  if (is.null(names)) names <- vars
  nvars  <- length(vars)
  if (!nvars) return(data)
  if (nvars != length(names)) {
    stop("ERROR with vars/names field in file.list")
  }
  typeFlag <- !is.null(type)
  if (typeFlag) {
    type <- toupper(type)
    if (length(type) == 1) type <- rep(type, times=nvars)
    if (nvars != length(type)) {
      stop("ERROR with vars/type field in file.list")
    }
  }
  nid <- length(id.var)
  if (length(id) != nid) {
    stop("ERROR in addColumn: length(id.var) != length(file.list$id.var)")
  }

  if (!is.data.frame(data)) {
    dfFlag <- 0
    if (typeFlag) {
      if (length(unique(type)) > 1) typeFlag <- 0
    }
  } else {
    dfFlag <- 1
  }

  # Read in the file
  if (file.list$file.type %in% c(3, 6, 8)) {
    file.list$method       <- 2
    file.list$what         <- "character"
    file.list$returnMatrix <- 1
    file.list$include.row1 <- file.list$header
  }
  if (is.data.frame(file.list$file)) {
    x <- file.list$file
    file.list$file <- NULL
    temp <- gc()
  } else {
    temp <- checkVars(file.list, file.list$vars)
    x <- loadData(file.list$file, file.list)
  }

  # Check id
  if (nid == 1) {
    if (length(unique(x[, id])) != nrow(x)) {
      print("WARNING: ids are not unique in file")
    }
    if (length(unique(data[, id.var])) != nrow(data)) {
      print("WARNING: ids are not unique in data")
    }
  } else {
    # Compare the id variables
    if (length(unique(data[, id.var[2]])) > length(unique(data[, id.var[1]]))) {
      # Switch
      temp      <- id.var[1]
      id.var[1] <- id.var[2]
      id.var[2] <- temp
      temp      <- id[1]
      id[1]     <- id[2]
      id[2]     <- temp
    } 
  }

  # Check variable names
  dimx <- dim(x)
  temp <- colnames(x)
  check.vec(vars, "file.list$vars", list(checkList=temp, maxValue=dimx[2], minValue=1))

  for (v in c(id, vars)) x[, v] <- unfactor(x[, v])
  for (v in id.var) data[, v] <- unfactor(data[, v])

  # Remove leading zeros
  if (op$leading0) {
    for (v in id.var) {
      data[, v] <- removeLeading0(data[, v])
    }
    for (v in id) {
      x[, v] <- removeLeading0(x[, v])
    }
  }

  # Get unique id.var[2]
  if (nid > 1) {
    levels <- unique(data[, id.var[2]])
  } else {
    levels <- 1
  }
  ret     <- NULL
  replace <- op$replace
  for (lev in levels) {
    if (nid > 1) {
      temp <- data[, id.var[2]] == lev
      d2   <- removeOrKeepRows(data, temp, which=1)
      temp <- x[, id[2]] == lev
      x2   <- removeOrKeepRows(x, temp, which=1)
    } else {
      d2 <- data
      x2 <- x
    } 

    # Match ids
    temp <- x2[, id[1]] %in% d2[, id.var[1]]
    x2   <- removeOrKeepRows(x2, temp, which=1)
    rows <- match(d2[, id.var[1]], x2[, id[1]])
    temp <- !is.na(rows)
    rows <- rows[temp]
    nr   <- length(rows)

    # For a matrix, define a matrix to combine
    if (!dfFlag) {
      i <- matrix(data=op$initValue, nrow=nrow(d2), ncol=nvars)
      colnames(i) <- names
      d2 <- cbind(d2, i)
    }

    # Add the columns
    for (i in 1:nvars) {
      if ((replace) && (dfFlag)) d2[, names[i]] <- op$initValue
      if (nr) {
        vec <- x2[rows, vars[i]] 
        if (typeFlag) {
          if (type[i] == "N") {
            vec <- as.numeric(vec)
          } else if (type[i] == "C") {
            vec <- as.character(vec)
          }
        }
        d2[temp, names[i]] <- vec
      }
    }
    ret <- rbind(ret, d2)

  } # END: for (lev in levels)

  ret

} # END: addColumn

# Function to order columns in a matrix or data frame
orderVars <- function(data, order) {

  # data     matrix or data frame with column names
  # order    Character vector

  if (ncol(data) == 1) return(data)
  cnames <- colnames(data)
  temp   <- order %in% cnames
  order  <- order[temp]
  temp   <- !(cnames %in% order)
  if (any(temp)) order <- c(order, cnames[temp])
    
  data <- data[, order]
  data

} # END: orderVars

# Function to remove leading zeros in character vectors
removeLeading0 <- function(vec) {

  # vec    Character vector

  if (is.numeric(vec)) return(vec)

  ret    <- vec
  numvec <- as.numeric(vec)
  temp   <- !is.na(numvec)
  if (any(temp)) ret[temp] <- as.character(numvec[temp])
  ret

} # END: removeLeading0


# Function to check that the specified delimiter at leasts exists in the first
#   row of the file.
checkDelimiter <- function(file.list) {

  file.list <- default.list(file.list, c("file", "file.type", "delimiter"),
                 list("ERROR", "ERROR", "ERROR"), error=c(1, 1, 1))

  if (file.list$file.type == 1) return(1)

  # Open file
  fid <- getFID(file.list$file, file.list)

  # Read 1 row
  x <- scan(fid, what="character", nlines=1, sep="\n")
  close(fid)

  ret <- grep(file.list$delimiter, x)
  if (!length(ret)) ret <- 0

  ret

} # END: checkDelimiter

# Function to remove leading/trailing white space
removeWhiteSpace <- function(str, leading=1, trailing=1) {

  if ((leading) && (trailing)) {
    ret <- gsub("^\\s+|\\s+$", "", str, perl=TRUE)
  } else if (leading) {
    ret <- gsub("^\\s+", "", str, perl=TRUE)
  } else if (trailing) {
    ret <- gsub("\\s+$", "", str, perl=TRUE)
  } else {
    ret <- str
  }

  ret

} # END: removeWhiteSpace

# Function for changing levels in a matrix or data frame
changeLevels.var <- function(data, var, old.levels, new.level,
                             new.var=NULL) {

  # data
  # var           Variable name or column number.
  #               Use NULL if data is a vector
  # old.levels    Current level(s) to change to new.level
  # new.level
  # new.var       New variable name or the old one if NULL
 
  # NA works with %in%

  if (is.null(dim(data))) {
    temp <- data %in% old.levels
    temp[is.na(temp)] <- FALSE
    data[temp] <- new.level
  } else {
    if (is.null(new.var)) new.var <- var
    temp <- data[, var] %in% old.levels
    temp[is.na(temp)] <- FALSE
    data[temp, new.var] <- new.level
  }

  data  

} # END: changeLevels.var

# Function for debugging 
debug.time <- function(time0, str=NULL) {

  if (!is.null(time0)) print(proc.time()-time0)
  if (!is.null(str)) print(str)
  return(proc.time())

} # END: debug.time

# Function to compute all interactions between 2 sets of vars
getInteractions <- function(data, vars1, vars2) {

  vars1 <- unique(vars1)
  vars2 <- unique(vars2)
  nv1   <- length(vars1)
  nv2   <- length(vars2)
  nv    <- nv1*nv2
  flip  <- 0
  if ((nv2 == 1) && (nv1 > 1)) {
    # Flip vars for efficiency
    temp  <- vars1
    vars1 <- vars2
    vars2 <- temp
    temp  <- nv1
    nv1   <- nv2
    nv2   <- temp
    flip  <- 1
  }
 
  new     <- matrix(data=NA, nrow=nrow(data), ncol=nv)
  newVars <- character(nv)
  if ((nv1 == 1) && (nv2 == 1)) {
    new[, 1] <- makeVector(data[, vars1])*makeVector(data[, vars2])
    newVars  <- paste(vars1, ".", vars2, sep="") 
  } else {
    start <- 1
    stop  <- nv2
    for (i in 1:nv1) {
       new[, start:stop] <- matrixMultVec(as.matrix(data[, vars2]), 
                              makeVector(data[, vars1[i]]), by=2)
       if (flip) {
         newVars[start:stop] <- paste(vars2, ".", vars1[i], sep="")
       } else {
         newVars[start:stop] <- paste(vars1[i], ".", vars2, sep="")
       }
       start <- stop + 1
       stop  <- stop + nv2
    }
  }

  colnames(new) <- newVars 
  data <- cbind(data, new) 
  
  list(data=data, newVars=newVars)

} # END: getInteractions

# Function to check for an error with try function
checkTryError <- function(obj, conv=1) {

  classObj <- class(obj)
  if ("try-error" %in% classObj) return(1)
  ret <- 0  

  # Check for convergence
  if (conv) {
    if (("glm" %in% classObj) || ("lm" %in% classObj)) {
      ret <- 1 - obj$converged
    } else if ("vglm" %in% classObj) {
      temp <- try(obj@criterion$loglikelihood, silent=TRUE)
      if ("try-error" %in% classObj) return(1)
      if ((!is.finite(temp)) || (is.null(temp))) ret <- 1
    } else if ("snp.logistic" %in% classObj) {
      if (is.null(obj$UML)) return(1)
    }
  }
  
  ret

} # END: checkTryError

# Function to check a list of type file.list
check.file.list <- function(flist, op=NULL) {

  # op          List with names 
  #  exist
  #  vars

  if (!is.list(flist)) flist <- list(file=flist)
  flist <- default.list(flist, c("file"), list("ERROR"), 
                        error=c(1))
  if (is.null(flist[["file.type", exact=TRUE]])) {
    flist$file.type <- getFileType(flist$file)
  }
  if (is.null(flist[["delimiter", exact=TRUE]])) {
    flist$delimiter <- getFileDelim(flist$file, type=flist$file.type)
  }
  if (is.null(flist[["header", exact=TRUE]])) {
    flist$header <- getFileHeader(flist)
  }
 
  op <- default.list(op, c("exist"), list(1)) 
  if (op$exist) {
    if (check.files(flist$file)) stop()
  }
  slist <- op[["subsetData", exact=TRUE]]
  if (!is.null(slist)) {
    slist <- check.subsetData.list(slist) 
    svars <- getSubsetDataVars(slist) 
  } else {
    svars <- NULL
  }
  vars <- op[["vars", exact=TRUE]]
  vars <- unique(c(vars, svars))
  if (!is.null(vars)) checkVars(flist, vars) 

  flist

} # END: check.file.list

# Function to return the variables in a list of type subsetData
getSubsetDataVars <- function(slist) {

  n   <- length(slist)
  ret <- character(n)
  for (i in 1:n) {
    ret[i] <- slist[[i]]$var
  }

  ret

} # END: getSubsetDataVars

# Function to check a list of type subsetData
check.subsetData.list <- function(slist) {

  n <- length(slist)
  for (i in 1:n) {
    temp <- default.list(slist[[i]], c("var", "operator", "value"), 
                         list("ERROR", "ERROR", "ERROR"), error=c(1, 1, 1))
  }

  slist

} # END: check.subsetData.list

# Function to normalize variable names
normVarNames <- function(cvec, op=NULL) {
 
  ret <- cvec
  ret <- gsub(">=", "_GTEQ_", ret, perl=TRUE)
  ret <- gsub("<=", "_LTEQ_", ret, perl=TRUE)
  ret <- gsub("==", "_EQ_", ret, perl=TRUE)
  ret <- gsub(">", "_GT_", ret, perl=TRUE)
  ret <- gsub("<", "_LT_", ret, perl=TRUE)
  ret <- gsub("%in%", "_IN_", ret, perl=TRUE)
  ret <- gsub("!=", "_NEQ_", ret, perl=TRUE)
  ret <- gsub(" ", ".", ret, perl=TRUE)

  str <- "[~`'!@#$%^&*()-+={}|\ ;<>?//]"
  ret <- gsub(str, "", ret, perl=TRUE)
  ret

} # END: normVarNames

# Function to get the SNP names from a matrix or data frame
getSnpNames <- function(obj, str="^rs[0-9]+$") {

  # obj    
  # str     PERL regular expression

  cnames <- colnames(obj)
  temp   <- grep(str, cnames, perl=TRUE)
  cnames <- cnames[temp]
  cnames

} # END: getSnpNames

# Function to get a unique variable name
getUniqueVarName <- function(vec, alen=4, nlen=4) {

  vec1 <- c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o",
            "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z")
  vec2 <- 0:9

  while (1) {
    t1  <- paste(sample(vec1, alen, replace=TRUE), collapse="", sep="")
    t2  <- paste(sample(vec2, alen, replace=TRUE), collapse="", sep="")
    ret <- paste(t1, t2, sep="")
    if (!(ret %in% vec)) return(ret)
  }

} # END: getUniqueVarName

