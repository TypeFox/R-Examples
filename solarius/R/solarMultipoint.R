#' Run multipoint linkage analysis.
#'
#' The linkage analysis is conducted in the following sequence:
#' parse input files of multipoint IBD matrices,
#' export data to a directory by \code{\link{df2solar}} function, 
#' run the polygenic analysis in a directory,
#' run the linkage analysis on the top of the polygenic analysis,
#' parse output files and 
#' store results in an object of \code{solarMultipoint} class (see \code{\link{solarMultipointClass}}).
#'
#' @note
#'  The user is responsible for computation of the MIBD matrices 
#'  required by \code{solarMultipoint} function via \code{mibddir} argument.
#'  Estimation of the IBD matrices is out of the scope of \code{solarius} package.
#'  More information is available in \url{http://solar.txbiomedgenetics.org/doc/05.chapter.html}.
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
#'@param mibddir 
#'    A character, the name of directory with files representing the IBD matrices in SOLAR format.
#'    These matrices can be evaluated by SOLAR or by other program.
#'    See SOLAR help for \code{ibd} command for more details
#'    (\url{http://solar.txbiomedgenetics.org/doc/91.appendix_1_text.html#ibd}).
#'@param chr 
#'    A character or one that can be coerced to that class,
#'    the value to be passed to SOLAR command \code{chromosome}.
#'    If it is missing, the default value is \code{"all"}.
#'    See SOLAR help for \code{chromosome} command for more details
#'    (\url{http://solar.txbiomedgenetics.org/doc/91.appendix_1_text.html#chromosome}).
#'@param interval
#'    A character or one that can be coerced to that class,
#'    the value to be passed to SOLAR command \code{interval}.
#'    If it is missing, the default value is \code{1}.
#'    See SOLAR help for \code{interval} command for more details.
#'@param multipoint.options
#'    A character of options to be passed to \code{multipoint} SOLAR command.
#'    For example, one or more LOD scores may be specified
#'    to run multi-pass linkage analysis.
#'    If the argument is \code{"3"}, that means a two-pass linkage scan
#'    with condition of the highest LOD of \code{3} to start the second pass.
#'    See SOLAR help page for \code{multipoint} command for more details
#'    (\url{http://solar.txbiomedgenetics.org/doc/91.appendix_1_text.html#multipoint}).
#'    The default value is \code{""}.
#'@param multipoint.settings
#'    A vector of characters, that contains SOLAR commands to be executed just before calling \code{multipoint}.
#'    For example, the fine mapping for the linkage scan can be disabled
#'    by setting the given argument to \code{"finemap off"}.
#'    The default value is \code{""}.
#'@param cores
#'    A positive integer, the number of cores for parallel computing.
#'    The default value is taken from \code{getOption("cores")}.
#'    If the default value is \code{NULL} then the number of cores is \code{1}.
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
#' ### load phenotype data 
#' data(dat30)
#'
#' ### load marker data
#' mibddir <- system.file('extdata', 'solarOutput',
#'   'solarMibdsCsv', package = 'solarius') 
#' list.files(mibddir)
#'
#' \dontrun{
#' ### basic (univariate) linkage model
#' mod <- solarMultipoint(trait1 ~ 1, dat30,
#'   mibddir = mibddir, chr = 5)
#' mod$lodf # table of results (LOD scores) (the highest 3.56)
#'
#' ### basic (bivariate) linkage model
#' mod <- solarMultipoint(trait1 + trait2 ~ 1, dat30,
#'   mibddir = mibddir, chr = 5)
#' mod$lodf # table of results (LOD scores) (the highest 2.74)
#'
#' ### two-pass linkage model
#' mod <- solarMultipoint(trait1 ~ 1, dat30,
#'   mibddir = mibddir, chr = 5,
#'   multipoint.options = "3")
#' mod$lodf # table of results (LOD scores, 1 pass) (the highest 2.74)
#' mod$lodf2 # table of results (LOD scores, 2 pass) (all nearly zero LOD scores)
#'
#' }
#' @export
solarMultipoint <- function(formula, data, dir,
  kinship,
  traits, covlist = "1",
  # input data to multipoint 
  mibddir,
  chr, interval, 
  multipoint.options = "", multipoint.settings = "",
  # output data from multipoint
  # misc
  cores = getOption("cores"),
  ...,
  verbose = 0) 
{
  tsolarMultipoint <- list()
  tsolarMultipoint$args <- proc.time()
  
  ### step 1: process par & create `out`
  mc <- match.call()
  
  # missing parameters
  stopifnot(!missing(mibddir))
  stopifnot(file.exists(mibddir))
  mibddir <- normalizePath(mibddir)
  
  # `chr`
  if(missing(chr)) {
    chr <- "all"
  }
  chr.str <- paste(chr, collapse = ",")
  
  # `interval`
  if(missing(interval)) {
    interval <- 1
  }

  
  # cores
  if(is.null(cores)) {  
    cores <- 1
  }
  
  is.tmpdir <- missing(dir)
  
  ### step 2.1: rename
  renames <- matchIdNames(names(data), skip.sex = TRUE) # match IDs, skip SEX
  unrelated <- ifelse(all(c("FA", "MO") %in% renames), FALSE, TRUE)
  
  # rename
  data <- rename(data, renames)
  
  ### step 2.2: take care of IDs in (1) pehnotypes in `data` argument; 
  # (2) IBD matrices in `mibddir` argument
  mibd.info <- get_info_mibd(mibddir)
  
  # match ids
  if(mibd.info$mibddir.format == "csv") {
    ids.data <- data$ID
    stopifnot(length(ids.data) > 0)
    ids.data <- as.character(ids.data)
    
    stopifnot(!is.null(mibd.info$mibddir.ids))
    ids.mibd <- mibd.info$mibddir.ids
    stopifnot(length(ids.mibd) > 0)
    
    # compare IDs
    ids.out <- ids.data[!(ids.data %in% ids.mibd)]
    if(length(ids.out) > 0) {
      warning(paste0(" Some individuals (IDs: ", 
        paste(ids.out, collapse = ", "), ") are not presented in IBDs. ",
        "Attempted to remove these individuals (FA/MO fields included) and pass new pedigrees to SOLAR."))

remove_ids_phen <- function(df, ids)
{
  stopifnot(all(ids %in% df$ID)) 
  
  # remove rows with `ID %in% ids`
  ind <- which(with(df, ID %in% ids))
  df <- df[-ind, ]
  
  # clean `FA` and `MO`
  # SOLAR ERROR: both parents must be known or unknown
  ids.fa <- with(df, as.character(FA) %in% ids)
  ids.mo <- with(df, as.character(MO) %in% ids)
  ids.par <- ids.fa | ids.mo

  df <- within(df, {
    FA[ids.par] <- "0"
    MO[ids.par] <- "0"  
  })
  
  return(df)
}

      data <- remove_ids_phen(data, ids.out)
    }
  }
  
  ### step 2.3: SOLAR dir
  if(is.tmpdir) {
    dir <- tempfile(pattern = "solarMultipoint-")
  }
  if(verbose > 1) cat(" * solarMultipoint: parameter `dir` is missing.\n")
  if(verbose > 1) cat("  -- temporary directory `", dir, "` used\n")

  ### step 3: compute a polygenic model by calling `solarPolygenic`
  tsolarMultipoint$polygenic <- proc.time()
  out <- solarPolygenic(formula, data, dir,
    kinship, traits, covlist, ..., verbose = verbose)

  # make a copy of `dir`
  files.dir <- list.files(dir, include.dirs = TRUE, full.names = TRUE)
  dir.poly <- file.path(dir, "solarPolygenic")
  stopifnot(dir.create(dir.poly, showWarnings = FALSE, recursive = TRUE))
  stopifnot(file.copy(from = files.dir, to = dir.poly, recursive = TRUE))
  
  ### step 4: add multipoint-specific slots to `out`
  tsolarMultipoint$premultipoint <- proc.time()
  
  out$multipoint <- list(call = mc,
    cores = cores,
    # input/output files
    dir.poly = dir.poly,
    out.dirs = "multipoint", out.chr = chr.str,
    # input/output data for multipoint
    mibddir = mibddir, mibd.info = mibd.info,
    chr = chr, chr.str = chr.str, interval = interval,
    multipoint.options = multipoint.options, multipoint.settings = multipoint.settings,
    tprofile = list(tproc = list()))

 ### step 5: prepare data for parallel computation (if necessary)
  out <- prepare_multipoint_files(out, dir)
  
  ### step 6: run multipoint
  tsolarMultipoint$runmultipoint <- proc.time()
  out <- run_multipoint(out, dir)
  
  ### clean 
  unlink(dir.poly, recursive = TRUE)
  
  if(is.tmpdir) {
    unlink(dir, recursive = TRUE)
    if(verbose > 1) cat("  -- solarMultipoint: temporary directory `", dir, "` unlinked\n")
  }

  ### return
  tsolarMultipoint$return <- proc.time()
  out$multipoint$tprofile$tproc$tsolarMultipoint <- tsolarMultipoint

  out$multipoint$tprofile <- try({
    procc_tprofile(out$multipoint$tprofile)
  }, silent = TRUE)
  
  oldClass(out) <- c("solarMultipoint", oldClass(out))
  return(out)
}

prepare_multipoint_files <- function(out, dir)
{  
  stopifnot(!is.null(out$multipoint$cores))
  stopifnot(!is.null(out$multipoint$out.dirs))
  stopifnot(!is.null(out$multipoint$out.chr))

  cores <- out$multipoint$cores
  out.dirs0 <- out$multipoint$out.dirs
  out.chr0 <- out$multipoint$out.chr
  chr <- out$multipoint$chr
  
  if(cores == 1) {
    out.dirs <- out.dirs0
    out.chr <- out.chr0
  } else {
    out.dirs <- paste(out.dirs0, 1:length(chr))
    out.chr <- chr
  }

  out$multipoint$out.dirs <- out.dirs
  out$multipoint$out.chr <- out.chr
  
  return(out)
}


run_multipoint <- function(out, dir)
{ 
  trun_multipoint <- list()
  trun_multipoint$args <- proc.time()
  
  cores <- out$multipoint$cores
  
  out.dirs <- out$multipoint$out.dirs
  out.chr <- out$multipoint$out.chr  
  
  ### run
  trun_multipoint$solar <- proc.time()
  parallel <- (cores > 1)
  if(parallel) {
 # load required R package doParallel
    stopifnot(requireNamespace("doParallel", quietly = TRUE))
    
    doParallel::registerDoParallel(cores = cores)
  }

  if(cores == 1) {
    out.gr <- llply(1, function(i) {
      solar_multipoint(dir, out, out.dirs, out.chr)
    })
  } else {
    num.out.dirs <- length(out.dirs)
    out.gr <- llply(1:num.out.dirs, function(i) {
      if(out$verbose) cat(" * solarMultipoint: ", i, "/", num.out.dirs, "batch...\n")
      solar_multipoint(dir, out, out.dirs[i], out.chr[i])
    }, .parallel = parallel)
  }
      
  ### process results
  trun_multipoint$results <- proc.time()
  results.list <- llply(out.gr, function(x) try({
    try(read_multipoint_lod(x$tab.dir, num.traits = length(out$traits)))
    #pedlod <- try(read_multipoint_pedlod(dir.multipoint))
  }))
  
  # -- try to union `lodf`
  num.passes <- try({
    results.list[[1]]$num.passes
  })
  
  lodf <- results.list 
  ret <- try({
    ldply(lodf, function(x) x$df)
  })

  if(class(ret)[1] != "try-error") {
    lodf <- ret
  }
  
  # passes 2 & 3 (if necessary)
  lodf2 <- data.frame()
  lodf3 <- data.frame()  
  if(class(num.passes) == "integer") {
    # pass 2
    if(num.passes >= 2) {
      lodf2 <- results.list
      ret <- try({
        ldply(lodf2, function(x) x$df2)
      })

      if(class(ret)[1] != "try-error") {
        lodf2 <- ret
      }
    }
    # pass 3
    if(num.passes >= 3) {
      lodf3 <- results.list
      ret <- try({
        ldply(lodf3, function(x) x$df3)
      })

      if(class(ret)[1] != "try-error") {
        lodf3 <- ret
      }
    }    
  }  
  # -- extract multipoint. solar outputs
  multipoint.solar <- llply(out.gr, function(x) x$solar)

  out.multipoint.solar <- list(cmd = llply(multipoint.solar, function(x) x$cmd), 
      solar.ok = llply(multipoint.solar, function(x) x$solar.ok))
    
  # -- final output
  out$lodf <- lodf
  out$lodf2 <- lodf2
  out$lodf3 <- lodf3    
  out$multipoint$num.passes <- num.passes
  out$multipoint$solar <- out.multipoint.solar
  
  ### return
  trun_multipoint$return <- proc.time()
  out$multipoint$tprofile$tproc$trun_multipoint <- trun_multipoint
      
  return(out)
}
