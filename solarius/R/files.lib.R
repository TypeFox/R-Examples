#----------------------------------
# Phenotype files
#----------------------------------

#' Read plain-text table files phen and ped
#'
#' The function reads the phenotype and pedigree data in SOLAR format,
#' i.e. plain-text tables.
#'
#' @param phen.file
#'    A character, path to phen file.
#' @param sep.phen
#'    A character, the field separator in phen file.
#'    The default value is \code{","}.
#' @param ped.file
#'    (optional) A character, path to ped file.
#' @param sep.ped
#'    A character, the field separator in ped file.
#'    The default value is \code{","}.
#' @param header
#'    Logical, indicating whether the file contains 
#'    the names of the variables as its first line.
#'    The default value is TRUE.
#' @param stringsAsFactors
#'    logical, indicating whether character vectors to be converted to factors.
#'    The default value is FALSE.
#' @param id.unique
#'    logical, indicating whether the IDs of individuals must be unique.
#'    The default value is TRUE.
#' @param sex.optional
#'    logical, indicating whether the SEX field must be presented.
#'    The default value is TRUE if \code{ped.file} is specified,
#'    and it is FALSE otherwise.
#' @return
#'    A data frame of phenotype data merged from two phen and ped files (merged by ID filed).
#'
#' @export
readPhen <- function(phen.file, sep.phen = ",",
  ped.file, sep.ped = ",", 
  header = TRUE, stringsAsFactors = FALSE,
  id.unique = TRUE, sex.optional)
{
  stopifnot(!missing(phen.file) | !missing(ped.file)) 
  
  if(missing(sex.optional)) {
    sex.optional <- ifelse(missing(ped.file), FALSE, TRUE)
  }
  
  stopifnot(header)
  stopifnot(!stringsAsFactors)
  
  stopifnot(id.unique)
  
  ### read `phen` file
  if(!missing(phen.file)) {
    sep <- sep.phen
    phen1 <- read.table(phen.file, nrows = 1,
      sep = sep, header = header, stringsAsFactors = stringsAsFactors)
    new.names <- matchIdNames(names(phen1), sex.optional = sex.optional)
    old.names <- names(new.names)
  
    ind <- which(names(phen1) %in% old.names)
    colClasses <- rep(as.character(NA), ncol(phen1))
    colClasses[ind] <- "character"
  
    phen <- read.table(phen.file, colClasses = colClasses,
      sep = sep, header = header, stringsAsFactors = stringsAsFactors)
  
    phen <- rename(phen, new.names)
  
    if(id.unique) {
      stopifnot(!any(duplicated(phen$ID)))
    }  
  }
  
  ### read `ped` if necessary
  if(!missing(ped.file)) {
    stopifnot(file.exists(ped.file))
  
    sep <- sep.ped  
    ped <- read.table(ped.file, colClasses = "character",
      sep = sep, header = header, stringsAsFactors = stringsAsFactors)
   
    new.names <- matchIdNames(names(ped))
    ped <- rename(ped, new.names)

    if(id.unique) {
      stopifnot(!any(duplicated(ped$ID)))
    }  
    
    if(!missing(phen.file)) {
      # merge `phen` & `ped`
      stopifnot("ID" %in% new.names)
      new.names2 <- new.names[!new.names %in% c("ID")]
    
      ind <- which(names(phen) %in% new.names2)
      if(length(ind)) {
        phen <- phen[, -ind]
      }
      
      dat <- base::merge(phen, ped, by = "ID", all = TRUE)      
    }
  }
  
  ### return
  if(missing(ped.file)) {
    return(phen)
  } else if(missing(phen.file)) {
    return(ped)
  } else {
    return(dat)
  }
}

#----------------------------------
# Association
#----------------------------------

check_assoc_files_exist <- function(genocov.files, snplists.files, snpmap.files)
{
  # genocov.files
  if(!missing(genocov.files)) {
    status <- laply(genocov.files, file.exists)
    if(!all(status)) {
      stop("Error in `check_assoc_files_exist`: ", sum(!status), " `genocov.files` do not exist.")
    }

    status <- laply(genocov.files, function(x) file.info(x)$isdir)
    if(any(status)) {
      stop("Error in `check_assoc_files_exist`: ", sum(status), " `genocov.files` are directories.")
    }
  }
      
  # snplists.files
  if(!missing(snplists.files)) {
    status <- laply(snplists.files, file.exists)
    if(!all(status)) {
      stop("Error in `check_assoc_files_exist`: ", sum(!status), " `snplists.files` do not exist.")
    }

    status <- laply(snplists.files, function(x) file.info(x)$isdir)
    if(any(status)) {
      stop("Error in `check_assoc_files_exist`: ", sum(status), " `snplists.files` are directories.")
    }
  }

  # snpmap.files
  if(!missing(snpmap.files)) {
    status <- laply(snpmap.files, file.exists)
    if(!all(status)) {
      stop("Error in `check_assoc_files_exist`: ", sum(!status), " `snpmap.files` do not exist.")
    }

    status <- laply(snpmap.files, function(x) file.info(x)$isdir)
    if(any(status)) {
      stop("Error in `check_assoc_files_exist`: ", sum(status), " `snpmap.files` are directories.")
    }
  }  
  return(invisible())  
}
  
#----------------------------------
# Plink files
#----------------------------------

read_plink_raw <- function(plink.raw)
{
  tab <- fread(plink.raw, select = 2, stringsAsFactors = FALSE, colClasses = "character")
  ids <- tab$IID

  snpcovdata <- fread(plink.raw, drop = 1:6, header = TRUE, colClasses = "numeric", na.strings = "NA")
  snpcovdata <- as.matrix(snpcovdata)
    
  rownames(snpcovdata) <- ids
  return(snpcovdata)
}

read_plink_map <- function(plink.map)
{
  snpmap <- fread(plink.map, header = FALSE)
  stopifnot(ncol(snpmap) == 4)
  
  setnames(snpmap, 1:4, c("chr", "snp", "cM", "bp")) 

  return(snpmap)
}

read_plink_ped <- function(plink.ped, plink.map)
{
  tab <- fread(plink.ped, select = 2, stringsAsFactors = FALSE, colClasses = "character")
  ids <- tab$V2

  map <- fread(plink.map, select = 2, header = FALSE, colClasses = "character")
  snps <- map$V2

  ped <- fread(plink.ped, drop = 1:6, header = FALSE, colClasses = "character", na.strings = "0", data.table = FALSE)
  stopifnot(nrow(ped) == length(ids))
  stopifnot(ncol(ped) == 2 * length(snps))
  
  # transform alleles to genotype in a loop for all individuals
  snpdata <- matrix(as.character(NA), nrow = length(ids), ncol = length(snps))
  for(i in 1:length(ids)) {
    alleleA <- ped[i, seq(1, by = 2, length.out = length(snps))]
    alleleB <- ped[i, seq(2, by = 2, length.out = length(snps))]
    
    ind.na <- which(is.na(alleleA) | is.na(alleleB))
            
    genotype <- paste(alleleA, alleleB, sep = "/")
    if(length(ind.na)) {
      genotype[ind.na] <- as.character(NA)
    }
    
    snpdata[i, ] <- genotype
  }
  
  # add colnames/rownames
  rownames(snpdata) <- ids
  colnames(snpdata) <- snps
  
  return(snpdata)
}

  
