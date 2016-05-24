#----------------------------------
# Utils
#----------------------------------

procc_tprofile <- function(tprofile)
{
  tprocf <- list()
  
  num.tproc <- length(tprofile$tproc)
  for(i in 1:num.tproc) {
    df <- ldply(tprofile$tproc[[i]])
    
    df <- rename(df, c(.id = "unit"))
    
    elapsed <- NULL # R CMD check: no visible binding    
    df <- within(df, {
      elapsed.diff <- c(elapsed[-1] - elapsed[-nrow(df)], elapsed[nrow(df)] - elapsed[1])
    })
    df <- within(df, {
      elapsed.prop <- elapsed.diff / elapsed.diff[nrow(df)]
    })
          
    tprocf[[i]] <- df
  }
  names(tprocf) <- names(tprofile$tproc)
  
  ### return
  tprofile$tprocf <- tprocf
  
  tprofile$cputime.sec <- tail(tprofile$tprocf[[num.tproc]]$elapsed.diff, 1)
  
  return(tprofile)
}

#----------------------------------
# Read/Write Files
#----------------------------------

#' Read SOLAR output files in a directory
#'
#' The function is used, for example, to read the files
#' in the output directory of the polygenic model.
#'
#' @param dir
#'    A character, path to the directory.
#' @return
#'    A list with file contents acquired by \code{readLines} function.
#'
#' @export
solarReadFiles <- function(dir)
{
  stopifnot(!missing(dir))
  stopifnot(file.exists(dir))

  filenames <- list.files(dir)    
  files <- list.files(dir, full.names = TRUE)

  stopifnot(length(files) > 0)

  out <- list()    
  for(i in 1:length(files)) {
    out[[filenames[i]]] <- readLines(files[i])
  }
  
  return(out)
}

#pedindex.out                                          
# 5 IBDID                 IBDID                       I
# 1 BLANK                 BLANK                       C
# 5 FATHER'S IBDID        FIBDID                      I
# 1 BLANK                 BLANK                       C
# 5 MOTHER'S IBDID        MIBDID                      I
# 1 BLANK                 BLANK                       C
# 1 SEX                   SEX                         I
# 1 BLANK                 BLANK                       C
# 3 MZTWIN                MZTWIN                      I
# 1 BLANK                 BLANK                       C
# 5 PEDIGREE NUMBER       PEDNO                       I
# 1 BLANK                 BLANK                       C
# 5 GENERATION NUMBER     GEN                         I
# 1 BLANK                 BLANK                       C
# 1 FAMILY ID             FAMID                       C
# 2 ID                    ID                          C
read_pedindex <- function(pedindex.out, ids.unique = TRUE)
{
  ### CDE
  pedindex.cde <- paste(tools::file_path_sans_ext(pedindex.out), "cde", sep = ".")

  cf <- read.fwf(pedindex.cde, skip = 1, widths = c(2, 22, 28, 2),
    stringsAsFactors = FALSE)
  names(cf) <- c("LEN", "FIELDNAME", "FIELD", "LETTER")
  
  FIELD <- NULL # R CMD check: no visible binding
  cf <- mutate(cf,
    FIELD = gsub(" ", "", FIELD))
  
  # names & widths for `pedindex.out` file
  ind <- which(cf$FIELD != "BLANK")
  pnames <- cf$FIELD
  pwidths <- cf$LEN

  pf <- read.fwf(pedindex.out, widths = pwidths, colClasses = "character")

  pf <- pf[, ind]
  names(pf) <- pnames[ind]
  
  for(i in 1:ncol(pf)) {
    pf[, i] <- as.character(pf[, i])
  }  
  
  for(i in 1:ncol(pf)) {
    pf[, i] <- gsub(" ", "", pf[, i])
  } 
  
  if(ids.unique) {
    stopifnot(!all(duplicated(pf$ID)))
  }
  
  return(pf)
}

read_phi2_gz <- function(phi2.gz)
{
  kf <- read.table(gzfile(phi2.gz), colClasses = c("character", "character", "numeric", "numeric"))
  names(kf) <- c("IBDID1", "IBDID2", "phi2", "delta7")

  IBDID1 <- IBDID2 <- NULL # R CMD check: no visible binding
  kf <- mutate(kf,
    IBDID1 = gsub(" ", "", IBDID1),
    IBDID2 = gsub(" ", "", IBDID2))
#> head(kf)
#  IBDID1 IBDID2      phi2  delta7
#1      1      1 0.1560492  0.2574
#2      1      1 1.0000000  1.0000
#3      2      2 1.0000000  1.0000
#4      3      3 1.0000000  1.0000

  # get rid of first two lines, which may be duplicated
  if(with(kf, IBDID1[1] == IBDID2[1] & IBDID1[2] == IBDID2[2])) {
    kf <- kf[-1, ]
  }
  
  return(kf)
}

read_map <- function(file, remove.prefix = FALSE)
{
  ### read `chr` in the first line
  chr.line <- readLines(file, n = 1)
  
  # case 1: patters like `c1.1.500`
  chr <- strsplit(chr.line, "\\.")[[1]][1] # now it is like `c1`
  
  chr <- gsub("chrom", "", chr)
  chr <- gsub("chr", "", chr)
  chr <- gsub("c", "", chr)
  
  chr <- as.integer(chr)

  # case 2: patters like `c.22.1`
  if(is.na(chr)) { 
    chr <- strsplit(chr.line, "\\.")[[1]][2] # now it is like `1`
  
    chr <- as.integer(chr)
  }
    
  stopifnot(!is.na(chr))
    
  ### read 2 columns: snp name & position in bp
  tab <- data.table(SNP = character(0), pos = integer(0), chr = character(0))
  
  ret <- suppressWarnings(try({
    tab <- fread(file, skip = 1, header = FALSE, sep = " ")
  
    stopifnot(ncol(tab) == 2)
    colnames(tab) <- c("SNP", "pos")
    
    # manually remove `snp_` prefix in `SNP` column
    if(remove.prefix) {
      SNP <- NULL # fix `no visible binding`
      tab[, SNP := gsub("^snp_", "", SNP)]
    }
          
    ### add `chr` column
    tab <- data.table(tab, chr = as.character(chr))
  }, silent = TRUE))
  
  return(tab)
}

read_map_files <- function(files, cores = 1, ...)
{
  parallel <- (cores > 1)
  if(parallel) {
    # load required R package doParallel
    stopifnot(requireNamespace("doParallel", quietly = TRUE))
    
    doParallel::registerDoParallel(cores = cores)
  }
  
  map.list <- llply(files, function(f) try({
    read_map(f, ...)}), .parallel = parallel)
  
  # -- try to union `snpf` slots in `snpf.list`
  ret <- try({
    rbindlist(map.list)
  })
  
  map <- data.table(SNP = character(0), pos = integer(0), chr = character(0))
  if(class(ret)[1] != "try-error") {
    map <- ret
    
    stopifnot(ncol(map) == 3)
    stopifnot(all(c("SNP", "pos", "chr") %in% colnames(map)))
  }
  
  return(map)
}
#----------------------------------
# Kinship functions
#----------------------------------

kmat2phi2 <- function(kmat, dir, kin2.gz = "kin2.gz")
{
  kf <- kmat2kf(kmat)

  kf2phi2(kf, dir, kin2.gz = kin2.gz)  
}

kf2phi2 <- function(kf, dir, kin2.gz = "kin2.gz")
{
  #stopifnot(requireNamespace("gdata", quietly = TRUE))
  
  pedindex.out <- file.path(dir, "pedindex.out")
  pf <- read_pedindex(pedindex.out)
  
  kf <- kf_match_pedindex(kf, pf)

  knames2 <- c("ID1", "ID2", "phi2")
  #knames2 <- c("IBDID1", "IBDID2", "phi2")  
  stopifnot(knames2 %in% names(kf))
  kf2 <- subset(kf, select = knames2)
  kf2 <- rename(kf2, c(ID1 = "id1", ID2 = "id2", phi2 = "matrix1"))
  #kf2 <- rename(kf2, c(IBDID1 = "id1", IBDID2 = "id2", phi2 = "matrix1"))
  
  knames <- c("IBDID1", "IBDID2", "phi2")
  stopifnot(knames %in% names(kf))
  kf <- subset(kf, select = knames)

  phi2.gz <- file.path(dir, kin2.gz)
  
  # @ http://helix.nih.gov/Documentation/solar-6.6.2-doc/08.chapter.html#phi2
  # - The coefficients should begin in the fourteenth character column, 
  #   or higher, counting the first character column as number one. 
  #   That is why `width = c(10, 10, 10)`.
  #ret <- gdata::write.fwf(kf, gzfile(phi2.gz),
  #  rownames = FALSE, colnames = FALSE,
  #  sep = " ", width = c(10, 10, 10))

  #kf$d7 <- 1.0
  #kf <- mutate(kf,
  #  phi2 = sprintf("%.7f", phi2),
  #  d7 = sprintf("%.7f", d7)
  #)    
  
  #ret <- gdata::write.fwf(kf, gzfile(phi2.gz),
  #  rownames = FALSE, colnames = FALSE, justify = "right",
  #  sep = " ", width = c(5, 5, 10, 10))
  
  ### CSV format
  ord <- with(kf2, order(as.integer(id1), as.integer(id2)))
  kf2 <- kf2[ord, ]
  
  ret <- write.table(kf2, gzfile(phi2.gz), quote = FALSE,
    row.names = FALSE, col.names = TRUE, sep = ",")
  
  return(invisible())
}


#----------------------------------
# Checkers
#----------------------------------

check_var_names <- function(traits, covlist, dnames) {
  ### traits
  stopifnot(all(traits %in% dnames))
  
  ### covlist
  covlist2 <- covlist
  # filter out term 1
  ind <-  grep("^1$", covlist2)
  if(length(ind) > 0) { 
    covlist2 <- covlist2[-ind] 
  }
  # filter out interaction terms with `*`
  ind <-  grep("\\*", covlist2)
  if(length(ind) > 0) { 
    covlist2 <- covlist2[-ind] 
  }
  # filter out interaction terms with `#`
  ind <-  grep("\\#", covlist2)
  if(length(ind) > 0) { 
    covlist2 <- covlist2[-ind] 
  }
  # filter out power terms with "^"
  ind <-  grep("\\^", covlist2)
  if(length(ind) > 0) { 
    covlist2 <- covlist2[-ind] 
  }
  
  # filter out covariates like `c("sex(trait1)", "age(trait2)")`
  covlist3 <- covlist2
  if(length(covlist2)) {
    ind <- grep("\\(", covlist2)
    if(length(ind)) {
      covlist3[ind] <- laply(strsplit(covlist3[ind], "\\("), function(x) x[1])
    }
  }
  stopifnot(all(covlist3 %in% dnames))
  
  return(invisible())
}

#' Match ID column names
#'
#' The function automates the process of renaming ID fields in data frame of phenotypes
#' into a SOLAR format.
#'
#' @param fields
#'    The name of fields or colulm names in data frame, which are candidates to be ID fields.
#' @param sex.optional
#'    Logical, indicating if through an error in the case SEX field is not found in input fields.
#' @param skip.sex
#'    Logical, indicating if the search for SEX filds is completely skipped.
#' @return
#'    A named character vector, that can be directly passed to \code{rename} function of \code{plyr} package.
#'
#' @example inst/examples/example-fields.R
#' @export
matchIdNames <- function(fields, sex.optional = FALSE, skip.sex = FALSE)
{
  # `fields`: fields in data set
  # `names`: matched names
  
  find_name <- function(pat, fields) {
    names <- grep(pat, fields, value = TRUE)
    if(length(names) == 0) {
      stop("ID names not found; grep pattern '", pat, "'")
    }
    if(length(names) > 1) {
      stop("more than one ID names found (", paste(names, collapse = ", "), "); grep pattern '", pat, "'")
    }
    
    return(names)
  }
  
  ### option 1
  #out <- c("ID", "FAMID", "MO", "FA", "SEX")
  #out.names <- c(find_name("^id$|^ID$", fields),
  #  find_name("^famid$|^FAMID$", fields),
  #  find_name("^mo$|^MO$|^mother$|^MOTHER$", fields),
  #  find_name("^fa$|^FA$|^father$|^FATHER$", fields),
  #  find_name("^sex$|^SEX$", fields))
  #names(out) <- out.names

  ### option 2
  out <- c("ID")
  # ID field (obligatory)
  pat <- "^id$|^ID$"
  names <- grep(pat, fields, value = TRUE) 
  if(length(names) == 0) stop("ID name not found; grep pattern '", pat, "'")
  if(length(names) > 1)  stop("more than one ID names found (", paste(names, collapse = ", "), "); grep pattern '", pat, "'")
  out.names <- names
  # FAMID (optional)
  # -- pass 1
  pat <- "^famid$|^FAMID$"
  names <- grep(pat, fields, value = TRUE) 
  if(length(names) > 1)  stop("more than one FAMID names found (", paste(names, collapse = ", "), "); grep pattern '", pat, "'")   
  if(length(names) == 1) {
    out <- c(out, "FAMID")
    out.names <- c(out.names, names)
  } else {
    # -- pass 2
    pat2 <- "^fam$|^FAM$"
    names2 <- grep(pat2, fields, value = TRUE) 
    if(length(names2) > 1)  stop("more than one FAMID names found (", paste(names2, collapse = ", "), "); grep pattern '", pat2, "'")
    if(length(names2) == 1) {
      out <- c(out, "FAMID")
      out.names <- c(out.names, names2)
    }  
  }
  # MO (optional)
  pat <- "^mo$|^MO$|^mother$|^MOTHER$"
  names <- grep(pat, fields, value = TRUE) 
  if(length(names) > 1)  stop("more than one MO names found (", paste(names, collapse = ", "), "); grep pattern '", pat, "'")   
  if(length(names) == 1) {
    out <- c(out, "MO")
    out.names <- c(out.names, names)
  }
  # FA (optional)
  pat <- "^fa$|^FA$|^father$|^FATHER$"
  names <- grep(pat, fields, value = TRUE) 
  if(length(names) > 1)  stop("more than one FA names found (", paste(names, collapse = ", "), "); grep pattern '", pat, "'")   
  if(length(names) == 1) {
    out <- c(out, "FA")
    out.names <- c(out.names, names)
  }
  if(!skip.sex) {
    # SEX field (obligatory / optional if `sex.optional` argumnet is TRUE)
    pat <- "^sex$|^SEX$"
    names <- grep(pat, fields, value = TRUE) 
    if(!sex.optional) {
      if(length(names) == 0) stop("SEX name not found; grep pattern '", pat, "'")
    }
    if(length(names) > 1)  stop("more than one SEX names found (", paste(names, collapse = ", "), "); grep pattern '", pat, "'")
    if(length(names) == 1) {
      out <- c(out, "SEX")
      out.names <- c(out.names, names)
    }
  }
  # PROBND (optional)
  pat <- "^PROBND$|^probnd$"
  names <- grep(pat, fields, value = TRUE) 
  if(length(names) > 1)  stop("more than one PROBND names found (", paste(names, collapse = ", "), "); grep pattern '", pat, "'")   
  if(length(names) == 1) {
    out <- c(out, "PROBND")
    out.names <- c(out.names, names)
  }
  # HHID  (optional)
  pat <- "^HHID$|^hhid$"
  names <- grep(pat, fields, value = TRUE) 
  if(length(names) > 1)  stop("more than one HHID names found (", paste(names, collapse = ", "), "); grep pattern '", pat, "'")   
  if(length(names) == 1) {
    out <- c(out, "HHID")
    out.names <- c(out.names, names)
  }
  # MZTWIN (optional)
  pat <- "^MZTWIN$|^mztwin$"
  names <- grep(pat, fields, value = TRUE) 
  if(length(names) > 1)  stop("more than one MZTWIN names found (", paste(names, collapse = ", "), "); grep pattern '", pat, "'")   
  if(length(names) == 1) {
    out <- c(out, "MZTWIN")
    out.names <- c(out.names, names)
  }
  
  names(out) <- out.names    
  
  return(out)
}

#' Match map column names
#'
#' The function searches for fields that correspond to those 
#' in SOLAR format of marker map file.
#'
#' @param fields
#'    The name of fields or colulm names in data frame, which are candidates to be map fields.
#' @return
#'    A named character vector, that can be directly passed to \code{rename} function of \code{plyr} package.
#'
#' @examples
#'  data(dat50)
#'  matchMapNames(names(snpdata))
#' @export
matchMapNames <- function(fields)
{
  # `fields`: fields in data set
  # `names`: matched names
  
  ### matching
  out <- c("SNP")
  # ID field (obligatory)
  pat <- "^name$|^snp$|^SNP$|^Snp$|^marker$"
  names <- grep(pat, fields, value = TRUE) 
  if(length(names) == 0) stop("SNP name not found; grep pattern '", pat, "'")
  if(length(names) > 1)  stop("more than one SNP names found (", paste(names, collapse = ", "), "); grep pattern '", pat, "'")
  out.names <- names
  # pos 
  pat <- "^pos$|^position$|^Position$|^bp$"
  names <- grep(pat, fields, value = TRUE) 
  if(length(names) > 1)  stop("more than one `pos` names found (", paste(names, collapse = ", "), "); grep pattern '", pat, "'")   
  if(length(names) == 1) {
    out <- c(out, "pos")
    out.names <- c(out.names, names)
  }
  # chr 
  pat <- "^chr$|^chrom$|^chromosome$"
  names <- grep(pat, fields, value = TRUE) 
  if(length(names) > 1)  stop("more than one `chr` names found (", paste(names, collapse = ", "), "); grep pattern '", pat, "'")   
  if(length(names) == 1) {
    out <- c(out, "chr")
    out.names <- c(out.names, names)
  }

  names(out) <- out.names    
  
  return(out)
}
