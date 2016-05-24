
#----------------------------------
# Association functions
#----------------------------------

# function parameters
# - snplists.files
# - out.dir
# - out.file
#
# slots in the object `out` required to run SOLAR> mga
# - out$traits
# - out$solar$model.filename
# - out$solar$phe.filename)
solar_assoc <- function(dir, out, genocov.files, snplists.files, out.dir, out.file, tmp.dir = FALSE)
{
  ### check arguments
  stopifnot(file.exists(dir))
  
  ### var
  if(!tmp.dir) {
    dir.assoc <- dir
    out.path <- out.dir
    tab.file <- file.path(dir, out.path, out.file)
  } else {
    dir.poly <- out$assoc$dir.poly
    out.path <- file.path(dir, out.dir)

    snplists.files.local <- out$assoc$snplists.files.local
    genocov.files.local <- out$assoc$genocov.files.local
  
    if(snplists.files.local) {
      snplists.files <- file.path(dir, snplists.files)
    }
    if(genocov.files.local) {
      genocov.files <- file.path(dir, genocov.files)
    }

    ### check files exist  
    stopifnot(all(file.exists(genocov.files)))
    stopifnot(all(file.exists(snplists.files)))
    stopifnot(dir.create(out.path))

    ### normalize path
    genocov.files <- normalizePath(genocov.files)
    snplists.files <- normalizePath(snplists.files)
    out.path <- normalizePath(out.path)
  
    ### create a temp. dir.
    dir.assoc <- file.path(dir, paste0("solar_assoc_", basename(out.file)))
    stopifnot(dir.create(dir.assoc))
    files.dir <- list.files(dir.poly, include.dirs = TRUE, full.names = TRUE)
    stopifnot(file.copy(from = files.dir, to = dir.assoc, recursive = TRUE))
    
    tab.file <- file.path(out.path, out.file)
  }
  
  ### SOLAR options/settings
  assoc.options <- out$assoc$assoc.options
      
  ### make `cmd`
  # local variables
  trait.dir <- paste(out$traits, collapse = ".")
  model.path <- file.path(trait.dir, out$solar$model.filename)
  
  # combine variables into a command string `cmd`
  cmd <- c(
    # needed to reload phenotypes, in order to avoid errors like 
    # `Phenotype named snp_s1 found in several files`
    #  - diagnostic: go to SOLAR dir and show phenos by `pheno` command,
    #    and you will see `snp.genocov` duplicated
    #paste("load pheno", out$solar$phe.filename), 
    paste("load model", model.path),
    paste("outdir", out.path),
    # mga option `-files snp.genocov` is not passed, as that provokes pheno-dulicates 
    # (SOLAR's strange things)
    paste("mga ", "-files ", paste(genocov.files, collapse = " "), 
      " -snplists ", paste(snplists.files, collapse = " "), " -out ", out.file, 
      " ", paste(assoc.options, collapse = " "), sep = ""))

  ### run solar    
  ret <- solar(cmd, dir.assoc, result = FALSE) 
  # `result = FALSE`, because all assoc. results are printed to output
  
  solar.ok <- file.exists(tab.file)

  ### clean
  if(tmp.dir) {
    unlink(dir.assoc, recursive = TRUE)
  }
    
  ### return  
  out <-  list(solar = list(cmd = cmd, solar.ok = solar.ok), tab.file = tab.file)

  return(out)
}

#----------------------------------
# Annotation functions
#----------------------------------

#' Annotate SNPs 
#'
#' The function annotates SNPs based on \code{NCBI2R} R package,
#' in particular \code{AnnotateSNPList} function.
#'
#' See \url{https://ncbi2r.wordpress.com/} for more details.
#'
#' @param x 
#'    An object of class \code{solarAssoc} or a character vector of SNPs.
#' @param mode
#'    A character with the mode of SNPs selection.
#'    Possible values are \code{"significant"}, \code{"top"} and \code{"all"}.
#'    The default value is \code{"significant"}.
#' @param alpha
#'    A numeric value from 0 to 1, the significance level after Bonferroni multiple-test correction.
#'    Corresponds to \code{mode} equal to \code{"significant"}.
#' @param num.top
#'    An integer value, the number of top SNPs to be annotated.
#'    Corresponds to \code{mode} equal to \code{"top"}.
#'    The default value is 10.
#' @param query.size
#'    An integer, the maximum number of SNPs allowed for a single query to the NCBI database.
#'    The default value is \code{500}.
#'    See also the help page for \code{NCBI_snp_query} function in \code{rsnps} package.
#'    If the number of SNPs is greater than \code{query.size}, 
#'    then the query is split into batches automatically.
#' @param verbose
#'    An integer, the verbose level.
#'    The default value is \code{0}.
#' @param ...
#'    Additional arguments.
#' @return
#'    A data table with annotation results.
#'
#' @export
annotateSNPs <- function(x, mode = c("significant", "top", "all"), 
  alpha = 0.05,
  num.top = 10, 
  query.size = 500,
  verbose = 0)
{
  ### inc
  stopifnot(requireNamespace("rsnps"))
  
  ### args
  mode <- match.arg(mode)

  ### inc
  stopifnot(requireNamespace("rsnps", quietly = TRUE))
  
  ### get list of SNPs
  if(class(x)[1] == "solarAssoc") {
    num.snps <- nrow(x$snpf)
    if(mode == "significant") { 
      SNP <- pSNP <- NULL # # due to R CMD check: no visible binding 
      snplist <- subset(x$snpf, pSNP <= (alpha/num.snps))[, SNP]
    
      if(length(snplist) == 0) {  
        warning("There are no significant SNPs.")
        return(invisible())
      }
    } else if(mode == "top") {
      num.top <- min(num.top, num.snps)
      
      ord <- order(x$snpf$pSNP)
      ord <- ord[seq(1, num.top)]
      
      snplist <- x$snpf[ord, SNP]
    } else if(mode == "all") {
      snplist <- x$snpf[, SNP]
    } else {
      stop("`mode` is unknown.")
    }
    
    snplist <- as.character(snplist)
  } else if (class(x)[1] == "character") {
    snplist <- x
  } else {
     stop("Class of `x` is unknown.")
  }  
  
  ### annotate
  # fitler out not-rs SNPs
  snplist <- grep("^rs", snplist, value = TRUE)
  
  # break down `snplist` into batches of size 500, if it necessary
  if(length(snplist) > query.size) {
    out <- list()
    
    ind <- seq(1, length(snplist)) # [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20
    num.queries <- ceiling(length(snplist) / query.size) 
    jnd <- cut(ind, num.queries, labels = FALSE) #  [1] 1 1 1 1 1 1 1 2 2 2 2 2 2 3 3 3 3 3 3 3
    for(i in seq(1, num.queries)) {
      if(verbose) cat(" * query", i, "/", num.queries, "\n")
  
      out[[i]] <- try(rsnps::NCBI_snp_query(snplist[ind[jnd == i]]))
    }
    knd <- (laply(out, class) == "data.frame")
    out <- out[knd]

    annot <- ldply(out, rbind)
  } else {
    annot <- rsnps::NCBI_snp_query(snplist)
  }
  # > annot  
  #      Query Chromosome    Marker Class Gene Alleles Major Minor   MAF        BP
  #1 rs2731672          5 rs2731672   snp  F12     A/G     A     G 0.484 177415472
  
  return(annot)
}

annotateSNPs.NCBI2R <- function(x, mode = c("significant", "top", "all"), 
  alpha = 0.05,
  num.top = 10, 
  capture.output = FALSE,
  ...)
{
  mode <- match.arg(mode)

  fun.exists <- exists("AnnotateSNPList", mode = "function")
  if(!fun.exists) {
    warning("R package `NCBI2R` is required. Load the package before calling `annotateSNPs` function. See the web page https://ncbi2r.wordpress.com/ for more details.")
    return(invisible())
  }
  
  ### get list of SNPs
  if(class(x)[1] == "solarAssoc") {
    num.snps <- nrow(x$snpf)
    if(mode == "significant") { 
      SNP <- pSNP <- NULL # # due to R CMD check: no visible binding 
      snplist <- subset(x$snpf, pSNP <= (alpha/num.snps))[, SNP]
    
      if(length(snplist) == 0) {  
        warning("There are no significant SNPs.")
        return(invisible())
      }
    } else if(mode == "top") {
      num.top <- min(num.top, num.snps)
      
      ord <- order(x$snpf$pSNP)
      ord <- ord[seq(1, num.top)]
      
      snplist <- x$snpf[ord, SNP]
    } else if(mode == "all") {
      snplist <- x$snpf[, SNP]
    } else {
      stop("`mode` is unknown.")
    }
    
    snplist <- as.character(snplist)
  } else if (class(x)[1] == "character") {
    snplist <- x
  } else {
     stop("Class of `x` is unknown.")
  }  
  
  ### annotate
  cmd.annotate <- "AnnotateSNPList(snplist)"
  if(capture.output) {
    tmp.file <- tempfile("AnnotateSNPList-")
    capture.output({
      tab <- try({ 
        eval(parse(text = cmd.annotate))
      })}, 
      file = tmp.file)
    ret <- unlink(tmp.file)
  } else {
    tab <- try({ 
      eval(parse(text = cmd.annotate))
    }) 
  }
      
  if(class(tab) == "try-error") {
    warning("`AnnotateSNPList` function in R package `NCBI2R` has not been completed.")
    return(invisible())
  }

  tab[1:12]
}


#----------------------------------
# Support functions
#----------------------------------

guess_snpformat <- function(snpdata, snpfile.gen, snpfile.genocov)
{
  if(!missing(snpdata)) {
    # first row
    vals <- snpdata[1, ]
    vals.class <- class(vals)
  } else {
    stop("error in `guess_snpformat`: snp values not extracted.")
  }
    
  # - option 1: `/` format
  if(vals.class == "character") {
    format.ok <- all(grepl("\\/", vals))
    if(!format.ok) {
      stop("error in `guess_snpformat`: snp data (1st row) is character, but not in supported `/` format.")
    } else {
      return("/")
    }
  } else if(vals.class == "numeric") {
    # - option 2: `012` format
    format.ok <- all(vals %in% c(0, 1, 2))
    if(format.ok) {
      return("012")
    } 
    else {
      # - option 3: `0.1` format
      return("0.1")
    }
  }
}
