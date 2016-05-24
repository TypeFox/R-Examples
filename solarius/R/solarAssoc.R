#' Run association analysis.
#'
#' The association analysis is conducted in the following sequence:
#' parse input files of SNP markers,
#' export data to a directory by \code{\link{df2solar}} function, 
#' run the polygenic analysis in a directory,
#' run the association analysis on the top of the polygenic analysis,
#' parse output files and 
#' store results in an object of \code{solarAssoc} class (see \code{\link{solarAssocClass}}).
#'
#' @note
#'  \code{solarAssoc} function accepts input genetic data in three formats:
#'  SOLAR (\code{genocov.files}, \code{snplists.files}, \code{snpmap.files} and \code{param mga.files} arguments),
#'  R data frame or matrix (\code{snpdata}, \code{snpcovdata} and \code{snpmap} arguments),
#'  and plink (\code{plink.ped}, \code{plink.map} and \code{plink.raw} arguments).
#'
#'  For large-size problems, the user is recommended to prepare the genetic data in SOLAR format
#'  and to split them into batches of size, for example, 1,000 markers.
#'  The use of the other two R and plink formats is not optimized for large-scale scenarios.
#'
#'@param formula
#'  an object of class \code{formula} or one that can be coerced to that class.
#'  It is a symbolic description of fixed effects (covariates) to be fitted. 
#'@param data
#'    A data frame containing the variables in the model, 
#'    including ID fields needed to construct random effects: genetic and house-hold (both optional).
#'    Other classes such as list, environment or object coercible by \code{as.data.frame} to a data frame
#'    are not supported.
#'@param dir
#'    an optional character string, the name of directory,
#'    where SOLAR performs the analysis.
#'    In this case, the analysis within related input/output files is 
#'    conducted in the given folder instead of a temporary one
#'    (the default work flow).
#'@param kinship
#'    A matrix of the kinship coefficients (custom kinship matrix). 
#'    The IDs are required to be in row and column names.
#'@param traits
#'    a vector of characters to specify trait(s) in the model. It is alternative to the formula interface.
#'@param covlist 
#'    a vector of characters to specify fixed effects (covariates) in the model.
#'    It is alternative to the formula interface.
#'    The default value is \code{"1"}.
#'@param snpformat
#'    a character, the format of SNP data passsed by \code{snpdata} argument.
#'    Currently, this argument is not used.
#'@param snpdata
#'    A matrix of SNP data.
#'    SNPs are given in the columns, and individuals correspond to the rows. 
#'    The IDs of individuals are required to be in row and column names.
#'@param snpcovdata
#'    A matrix of SNP data, which are converted to covariates (numeric format).
#'    SNPs are given in the columns, and individuals correspond to the rows. 
#'    The IDs of individuals are required to be in row and column names.
#'@param snpmap
#'    A data.frame of annotation for SNPs.
#'@param snplist
#'    a vector of characters, the names of SNPs to be used in the analysis.
#'    This argument may be used when a subset of SNPs is of the interest.
#'@param snpind
#'    a vector of positive integers, the indices of SNPs to be used in the analysis.
#'    This argument may be used when a subset of SNPs is of the interest.
#'@param genocov.files
#'    A vector of characters, the file paths to \code{genocov} SOLAR files.
#'@param snplists.files
#'    A vector of characters, the file paths to \code{snplists} SOLAR files.
#'@param snpmap.files
#'    A vector of characters (optional), the file paths to \code{snpmap} SOLAR files.
#'@param mga.files
#'    A list with 2-3 elements, where each element is a vector of characters.
#'    This argument is an alternative to the other three
#'    \code{genocov.files}, \code{snplists.files} and \code{snpmap.files}.
#'    The element 3 of the list is optional.
#'@param plink.ped
#'    A character, the file path to genotype data in plink .ped format.
#'    Two columns are used per genotype.
#'@param plink.map
#'    A character, the file path to genotype annotation data in plink .map format.
#'@param plink.raw
#'    A character, the file path to genotype data in allele-dosage plink format 
#'    (an example plink command: \code{plink --noweb --file dat50  --recodeA}).
#'    One column is used per genotype. 
#'@param assoc.outformat
#'    A character, the output format.
#'    Possible values are \code{"df"}, \code{"outfile"} and \code{"outfile.gz"}.
#'    Currently, the only supported output format is \code{"df"}.
#'    That means the table of results is stored in \code{snpf} slot of a returned object.
#'@param assoc.outdir
#'    a character, the path to the output directory.
#'    Currently, this argument is not used.
#'@param assoc.options
#'    A character, specific options to be passed to \code{mga} SOLAR command.
#'@param cores
#'    A positive integer, the number of cores for parallel computing.
#'    The default value is taken from \code{getOption("cores")}.
#'    If the default value is \code{NULL} then the number of cores is \code{1}.
#'@param batch.size
#'    An integer, the number of SNPs per batch for parallel computation.
#'    The default value is \code{1000}.
#'@param ...
#'    Arguments to be passed to  \code{\link{solarPolygenic}} function.
#'    For example, one of such arguments may be \code{polygenic.settings = "option EnableDiscrete 0"}.
#'    Arguments of \code{solarMultipoint},
#'    which are also passed to \code{\link{solarPolygenic}},
#'    include \code{formula}, \code{data}, \code{dir}, 
#'    \code{kinship}, \code{traits} and \code{covlist}.
#'@param verbose 
#'    An non-negative integer of the verbose level.
#'    The default value is \code{0}.
#'
#'@examples
#' ### load data 
#' data(dat50)
#' dim(phenodata)
#' dim(kin)
#' dim(genodata)
#'
#' \dontrun{
#' ### basic (univariate) association model with a custom kinship
#' mod <- solarAssoc(trait~age+sex, phenodata,
#'   kinship = kin, snpdata = genodata)
#' mod$snpf # table of results for 50 SNPs
#' }
#' @export
solarAssoc <- function(formula, data, dir,
  kinship,
  traits, covlist = "1",
  # input data to association 
  snpformat, snpdata, snpcovdata, snpmap,
  snplist, snpind,
  genocov.files, snplists.files, snpmap.files,
  mga.files,
  plink.ped, plink.map, plink.raw, 
  # output data from association
  assoc.outformat = c("df", "outfile", "outfile.gz"), assoc.outdir, 
  # SOLAR options/settings
  assoc.options = "",
  # misc
  cores = getOption("cores"), batch.size = 1000,
  ...,
  verbose = 0) 
{
  tsolarAssoc <- list()
  tsolarAssoc$args <- proc.time()
  
  ### step 1: process par & create `out`
  mc <- match.call()
  
  ### check if input files exist
  ret <- check_assoc_files_exist(genocov.files, snplists.files, snpmap.files)
  
  # 1.1 missing parameters
  #if(missing(snpformat)) snpformat <- "012"
  
  missing.genocov.files <- missing(genocov.files)
  missing.snplists.files <- missing(snplists.files)
  missing.snpmap.files <- missing(snpmap.files)
  
  missing.mga.files <- missing(mga.files)

  missing.snpdata <- missing(snpdata)
  missing.snpmap <- missing(snpmap)
  missing.snpcovdata <- missing(snpcovdata)
  
  missing.plink.ped <- missing(plink.ped)
  missing.plink.map <- missing(plink.map)
  missing.plink.raw <- missing(plink.raw)

  # 1.2 process `mga.files` 
  if(!missing.mga.files & !missing.genocov.files) {
    stop("Error in `solarAssoc`: input SNP data must be given by either `mga.files` or `genocov.files` argument.")
  }
  
  if(!missing.mga.files) {
    stopifnot(class(mga.files) == "list")
    
    if(length(mga.files) < 2) {
      stop("Error in `solarAssoc`: `mga.files` must have 2 elements at least (genocov.files & snplists.files).")
    }
    # assign
    missing.genocov.files <- FALSE
    genocov.files <- mga.files[[1]]
    
    missing.snplists.files <- FALSE
    snplists.files <- mga.files[[2]]
    
    if(length(mga.files) > 2) {
      missing.snpmap.files <- FALSE
      snpmap.files <- mga.files[[3]]
    }
  }

  # 1.3 process `plink.raw`
  if(!missing.plink.raw & !missing.snpcovdata) {
    stop("Error in `solarAssoc`: input SNP covariate data must be given by either `plink.raw` or `snpcovdata` argument.")
  }
  
  if(!missing.plink.raw) {
    stopifnot(length(plink.raw) == 1)
    stopifnot(file.exists(plink.raw))
    
    # read data into `snpcovdata`
    missing.snpcovdata <- FALSE
    snpcovdata <- read_plink_raw(plink.raw)
  }

  # 1.4 process `plink.map`
  if(!missing.plink.map & !missing.snpmap) {
    stop("Error in `solarAssoc`: input SNP covariate data must be given by either `plink.map` or `snpmap` argument.")
  }
  
  if(!missing.plink.map) {
    stopifnot(length(plink.map) == 1)
    stopifnot(file.exists(plink.map))
    
    # read data into `snpcovdata`
    missing.snpmap <- FALSE
    snpmap <- read_plink_map(plink.map)
  }

  # 1.4 process `plink.ped`
  if(!missing.plink.ped & !missing.snpdata) {
    stop("Error in `solarAssoc`: input SNP covariate data must be given by either `plink.ped` or `snpdata` argument.")
  }
  
  if(!missing.plink.ped) {
    if(missing.plink.map) {
      stop("Error in `solarAssoc`: both `plink.ped` & `plink.map` must be specified.")
    }    
    stopifnot(length(plink.ped) == 1)
    stopifnot(file.exists(plink.ped))
    
    # read data into `snpcovdata`
    missing.snpdata <- FALSE
    snpdata <- read_plink_ped(plink.ped, plink.map)
  }
    
  # 1.5 check for input data argument
  if(missing.snpdata & missing.snpcovdata & missing.genocov.files) {
    stop("Error in `solarAssoc`: input SNP data must be given by either `snpdata`/`snpcovdata` or `genocov.files` argument.")
  }

  if(!missing.snpdata & !missing.snpcovdata) {
    stop("Error in `solarAssoc`: input SNP data must be given by either `snpdata` or `snpcovdata` or `genocov.files` argument.")
  }

  if(!missing.genocov.files) {
    if(missing.snplists.files) {
      stop("Error in `solarAssoc`: both genocov.files & snplists.files must be specified.")
    }
    if(length(genocov.files) > 1) {
      if(length(genocov.files) != length(snplists.files)) {
        stop("Error in `solarAssoc`: if several `genocov.files` (", length(genocov.files), 
          ") are given, then the same number of `snplists.files` (",length(snplists.files), ") is required.\n",
          "  --  SOLAR processes these files differently: `genocov.files` are processed one by one in the analysis chain, ",
          "while `snplists.files` are read all together before the analysis starts.\n",
          "  --  When several groups of files are given, SOLAR is called for each group independelty. ",
          "Groping, e.g. SNPS by chromosome, supposed to be under the user control (before SOLAR is called).")
      }
    }
  }

  assoc.informat <- ifelse(!missing.genocov.files, "genocov.file",
    ifelse(!missing.snpdata, "snpdata",
    ifelse(!missing.snpcovdata, "snpcovdata",
    stop("ifelse error in processing `assoc.informat`"))))
  if(assoc.informat == "genocov.file") {
    if(length(genocov.files) > 1) {
      assoc.informat <- "genocov.files"
    }
  }

  # check map files
  if(!missing.snpmap & !missing.snpmap.files) {
    stop("Error in `solarAssoc`: input SNP maps must be given by either `snpmap` or `snpmap.files` arguments.")
  }
  if(!missing.genocov.files & !missing.snpmap.files) {
    if(length(genocov.files) > 1) {
      if(length(genocov.files) != length(snpmap.files)) {
        stop("Error in `solarAssoc`: if several `genocov.files` (", length(genocov.files), 
          ") are given, then the same number of `snpmap.files` (",length(snplists.files), ") is required.")
      }
    }
  }  
  
  assoc.mapformat <- ifelse(!missing.snpmap.files, "snpmap.file",
    ifelse(!missing.snpmap, "snpmap", "default"))
  if(!missing.genocov.files & !missing.snpmap.files) {
    if(length(genocov.files) > 1) {
      if(assoc.mapformat == "snpmap.file") {
        stopifnot(length(genocov.files) == length(snpmap.files))
        assoc.mapformat <- "snpmap.files"
      }
    }
  }
    
  # check for matrix format  
  if(!missing.snpdata) {
    stopifnot(class(snpdata) == "matrix")
  }
  if(!missing.snpcovdata) {
    stopifnot(class(snpcovdata) == "matrix")
  }
  
  # check `snplist` / `snpind` format
  stopifnot(any(missing(snplist), missing(snpind)))
  
  assoc.snplistformat <- ifelse(all(missing(snplist), missing(snpind)), "default",
    ifelse(!missing(snplist), "snplist",
    ifelse(!missing(snpind), "snpind",
    stop("ifelse error in processing `snplistformat`"))))
  if(assoc.informat == "genocov.files") {
    if(assoc.snplistformat != "default") {
      stop("Error in `solarAssoc`: `snplist`/`snpind` is not allowed for `genocov.files` input format.")
    }
  }

  if(missing(snplist)) {
    snplist <-  character()
  }
  if(missing(snpind)) {
    snpind <- integer()
  }
  
  # check for output data arguments
  assoc.outformat <- match.arg(assoc.outformat)
    
  # cores
  if(is.null(cores)) {  
    cores <- 1
  }
  
  is.tmpdir <- missing(dir)
  
  # gues format
  #if(missing(snpformat)) {
  #  snpformat <- guess_snpformat(snpdata, snpfile.gen, snpfile.genocov)
  #}

  ### step 2: SOLAR dir
  if(is.tmpdir) {
    dir <- tempfile(pattern = "solarAssoc-")
  }
  if(verbose > 1) cat(" * solarAssoc: parameter `dir` is missing.\n")
  if(verbose > 1) cat("  -- temporary directory `", dir, "` used\n")

  ### step 3: compute a polygenic model by calling `solarPolygenic`
  tsolarAssoc$polygenic <- proc.time()
  out <- solarPolygenic(formula, data, dir,
    kinship, traits, covlist, ..., verbose = verbose)

  # make a copy of `dir`
  files.dir <- list.files(dir, include.dirs = TRUE, full.names = TRUE)
  dir.poly <- file.path(dir, "solarPolygenic")
  stopifnot(dir.create(dir.poly, showWarnings = FALSE, recursive = TRUE))
  stopifnot(file.copy(from = files.dir, to = dir.poly, recursive = TRUE))
  
  ### step 3.1: add assoc.-specific slots to `out`
  tsolarAssoc$preassoc <- proc.time()
  if(missing.genocov.files) {
    genocov.files.local <- TRUE
    genocov.files <- "snp.genocov"
  } else {
    genocov.files.local <- FALSE
    genocov.files <- normalizePath(genocov.files)
  }
  
  if(missing.snplists.files) {
    snplists.files.local <- TRUE
    snplists.files <- "snp.geno-list"
  } else {
    snplists.files.local <- FALSE
    snplists.files <- normalizePath(snplists.files)
  }
  if(missing.snpmap.files) {
    snpmap.files <- character(0)
  } else {
    snpmap.files <- normalizePath(snpmap.files)
  }
  
  out$assoc <- list(call = mc, #snpformat = snpformat,
    cores = cores, batch.size = batch.size,
    # input par
    snplist = snplist, snpind = snpind,
    # input/output files
    dir.poly = dir.poly,
    genocov.files = genocov.files, genocov.files.local = genocov.files.local,
    genolist.file = "snp.geno-list", 
    snplists.files = snplists.files, snplists.files.local = snplists.files.local,
    snpmap.files = snpmap.files,
    out.dirs = "assoc", out.files = "assoc.out",
    # input/output data for association
    assoc.informat = assoc.informat,
    assoc.outformat = assoc.outformat,
    assoc.snplistformat = assoc.snplistformat,
    assoc.mapformat = assoc.mapformat,
    tprofile = list(tproc = list()),
    # SOLAR options/settings
    assoc.options = assoc.options)

  ### step 4: add genotype data to `dir`
  #snpdata <- format_snpdata(snpdata, snpformat)

  # maps (load previously to loading snp data)
  # -- SOLAR does not use this info. in assoc. analysis 
  #    neither output to the results file
  #if(!missing.snpmap) {
  #  ret <- snpmap2solar(snpmap, dir)
  #}
  
  if(out$assoc$assoc.informat == "snpdata") {
    ret <- snpdata2solar(snpdata, dir)
  } else if(out$assoc$assoc.informat == "snpcovdata") {
    ret <- snpcovdata2solar(snpcovdata, dir, out = out)
  }
  
  ### number of snps
  num.snps <- as.integer(NA)
  if(out$assoc$assoc.informat %in% c("snpdata", "snpcovdata")) {
    snps <- readLines(file.path(dir, out$assoc$genolist.file))
    num.snps <- length(snps)
  }

  out$assoc$num.snps <- num.snps
    
  ### step 6: prepare data for parallel computation (if necessary)
  out <- prepare_assoc_files(out, dir)
  
  ### step 7: run assoc  
  tsolarAssoc$runassoc <- proc.time()
  out <- run_assoc(out, dir)
  
  ### step 8: set keys
  tsolarAssoc$keyresults <- proc.time()
  
  SNP <- NULL # R CMD check: no visible binding
  if(class(out$snpf)[1] == "data.table") {
    setkey(out$snpf, SNP)
  }
  
  ### step 9: try to add mapping information
  ret <- suppressWarnings(try({
  if(out$assoc$assoc.mapformat == "snpmap") {
    # read map
    tsolarAssoc$map <- proc.time()
    snpmap <- as.data.table(snpmap)
    
    renames <- matchMapNames(names(snpmap))
    snpmap <- rename(snpmap, renames)
    
    snpmap <- subset(snpmap, select = renames)
    setkey(snpmap, SNP)
    
    # annotate 
    tsolarAssoc$annotate <- proc.time()
    #out$snpf <- data.table:::merge.data.table(out$snpf, snpmap, by = "SNP", all.x = TRUE)
    out$snpf <- merge(out$snpf, snpmap, by = "SNP", all.x = TRUE)
  } else if(out$assoc$assoc.mapformat %in% c("snpmap.file", "snpmap.files")) {
    # read map
    tsolarAssoc$map <- proc.time()
    snpmap <- read_map_files(out$assoc$snpmap.files, cores = out$assoc$cores, remove.prefix = TRUE)
    
    renames <- matchMapNames(names(snpmap))
    snpmap <- rename(snpmap, renames)
    
    snpmap <- subset(snpmap, select = renames)
    setkey(snpmap, SNP)
    
    # annotate 
    tsolarAssoc$annotate <- proc.time()
    #out$snpf <- data.table:::merge.data.table(out$snpf, snpmap, by = "SNP", all.x = TRUE)
    out$snpf <- merge(out$snpf, snpmap, by = "SNP", all.x = TRUE)
  }
  }, silent = TRUE))
  
  ### clean 
  unlink(dir.poly, recursive = TRUE)
  
  if(is.tmpdir) {
    unlink(dir, recursive = TRUE)
    if(verbose > 1) cat("  -- solarAssoc: temporary directory `", dir, "` unlinked\n")
  }

  ### return
  tsolarAssoc$return <- proc.time()
  out$assoc$tprofile$tproc$tsolarAssoc <- tsolarAssoc

  out$assoc$tprofile <- try({
    procc_tprofile(out$assoc$tprofile)
  }, silent = TRUE)
  
  oldClass(out) <- c("solarAssoc", oldClass(out))
  return(out)
}

prepare_assoc_files <- function(out, dir)
{  
  stopifnot(!is.null(out$assoc$cores))
  stopifnot(!is.null(out$assoc$batch.size))
  stopifnot(!is.null(out$assoc$genocov.files))
  stopifnot(!is.null(out$assoc$out.dirs))
  stopifnot(!is.null(out$assoc$out.files))

  assoc.snplistformat <- out$assoc$assoc.snplistformat
  assoc.informat <- out$assoc$assoc.informat
  
  cores <- out$assoc$cores
  batch.size <- out$assoc$batch.size
  
  genocov.files <- out$assoc$genocov.files
  snplists.files0 <- out$assoc$snplists.files
  genolist.file0 <- out$assoc$genolist.file
  snplists.files.local <- out$assoc$snplists.files.local
  
  out.dirs0 <- out$assoc$out.dirs
  out.files0 <- out$assoc$out.files

  snplists.files.local <- out$assoc$snplists.files.local

  ### case 1
  if(assoc.informat == "genocov.files") {
    stopifnot(length(genocov.files) == length(snplists.files0))
    
    snplists.files <- snplists.files0
    num.gr <- length(genocov.files)
    
    out.dirs <- rep(as.character(NA), num.gr)
    out.files <- rep(as.character(NA), num.gr)

    for(k in 1:num.gr) {
      out.dirs[k] <- paste(out.dirs0, k, sep = "")
      out.files[k] <- paste(out.files0, k, sep = "")
    }
  ### case 2
  } else {
    #### sub-case: `snplistformat` is `snplist`/`snpind`
    if(assoc.snplistformat == "snplist") {
      snplists.files <- genolist.file0
      
      writeLines(out$assoc$snplist, file.path(dir, snplists.files))

      snplists.files0 <- snplists.files
      snplists.files.local <- TRUE
    } else if(assoc.snplistformat == "snpind") {
      snplists.files <- genolist.file0
    
      if(snplists.files.local) {
        snps <- unlist(llply(file.path(dir, snplists.files0), function(x) readLines(x)))  
      } else {
        snps <- unlist(llply(snplists.files0, function(x) readLines(x)))
      }
      num.snps <- length(snps)
      stopifnot(all(out$assoc$snpind <= num.snps))
      
      snplist <- snps[out$assoc$snpind]
      writeLines(snplist, file.path(dir, snplists.files))

      snplists.files0 <- snplists.files
      snplists.files.local <- TRUE
    }
    
    if(cores == 1) {
      snplists.files <- snplists.files0
      out.dirs <- out.dirs0
      out.files <- out.files0
    } else {
    ### case 3
    # - use `genolist.file` to split lists of SNPs
      ### number of snps
      if(snplists.files.local) {
        snps <- unlist(llply(file.path(dir, snplists.files0), function(x) readLines(x)))  
      } else {
        snps <- unlist(llply(snplists.files0, function(x) readLines(x)))
      }
      num.snps <- length(snps)
      
      bsize <- min(batch.size, ceiling(num.snps / cores))
      
      num.gr <- ceiling(num.snps / bsize)
      
      # cut
      gr <- cut(1:num.snps, breaks = seq(1, num.snps, length.out = num.gr + 1), include.lowest = TRUE)

      snplists.files <- rep(as.character(NA), num.gr)
      out.dirs <- rep(as.character(NA), num.gr)
      out.files <- rep(as.character(NA), num.gr)

      for(k in 1:nlevels(gr)) {
        snplists.files.k <- paste(genolist.file0, k, sep = "")
        out.dirs.k <- paste(out.dirs0, k, sep = "")
        out.files.k <- paste(out.files0, k, sep = "")
        
        gr.k <- levels(gr)[k]
        snps.k <- snps[gr %in% gr.k]
    
        writeLines(snps.k, file.path(dir, snplists.files.k))
        
        snplists.files[k] <- snplists.files.k
        out.dirs[k] <- out.dirs.k
        out.files[k] <- out.files.k
      }
      snplists.files.local <- TRUE
    }
  }
  stopifnot(all(!is.na(snplists.files)))
  stopifnot(all(!is.na(out.dirs)))
  stopifnot(all(!is.na(out.files)))

  out$assoc$snplists.files <- snplists.files
  out$assoc$out.dirs <- out.dirs
  out$assoc$out.files <- out.files
  out$assoc$snplists.files.local <- snplists.files.local

  return(out)
}

run_assoc <- function(out, dir)
{ 
  trun_assoc <- list()
  trun_assoc$args <- proc.time()
  
  cores <- out$assoc$cores
    
  genocov.files <- out$assoc$genocov.files
  snplists.files <- out$assoc$snplists.files
  out.dirs <- out$assoc$out.dirs
  out.files <- out$assoc$out.files
  
  ### run
  trun_assoc$solar <- proc.time()
  parallel <- (cores > 1)
  if(parallel) {
    # load required R package doParallel
    stopifnot(requireNamespace("doParallel", quietly = TRUE))
    
    doParallel::registerDoParallel(cores = cores)
  }

  ### case 1
  if(length(genocov.files) > 1) {
    num.genocov.files <- length(genocov.files)
    out.gr <- llply(1:num.genocov.files, function(i) {
      if(out$verbose) cat(" * solarAssoc: ", i, "/", num.genocov.files, "genocov.files...\n")
      solar_assoc(dir, out, genocov.files[i], snplists.files[i], out.dirs[i], out.files[i], tmp.dir = parallel)
    }, .parallel = parallel)
  ### case 2 (length(genocov.files) == 1)
  } else if(cores == 1) {
    out.gr <- llply(1, function(i) {
      solar_assoc(dir, out, genocov.files, snplists.files, out.dirs, out.files)
    })
  ### case 2 (length(genocov.files) == 2)
  } else {
    num.snplists.files <- length(snplists.files)
    out.gr <- llply(1:num.snplists.files, function(i) {
      if(out$verbose) cat(" * solarAssoc: ", i, "/", num.snplists.files, "snplists.files...\n")
      solar_assoc(dir, out, genocov.files, snplists.files[i], out.dirs[i], out.files[i], tmp.dir = parallel)
    }, .parallel = parallel)
  }
      
  ### process results
  trun_assoc$results <- proc.time()
  snpf.list <- llply(out.gr, function(x) try({
    fread(x$tab.file)}), .parallel = parallel)
  
  # -- try to union `snpf` slots in `snpf.list`
  snpf <- snpf.list 
  ret <- try({
    rbindlist(snpf)
  })

  if(class(ret)[1] != "try-error") {
    snpf <- ret
    snpf <- rename(snpf, c("p(SNP)" = "pSNP"))
  }
    
  # -- extract assoc. solar outputs
  assoc.solar <- llply(out.gr, function(x) x$solar)
    out.assoc.solar <- list(cmd = llply(assoc.solar, function(x) x$cmd), 
      solar.ok = llply(assoc.solar, function(x) x$solar.ok))
    
  # -- final output
  out.assoc <- list(snpf = snpf, solar = out.assoc.solar)

  ### assign
  out$snpf <- out.assoc$snpf
  out$assoc$solar <- out.assoc$solar
  
  ### return
  trun_assoc$return <- proc.time()
  out$assoc$tprofile$tproc$trun_assoc <- trun_assoc
      
  return(out)
}
