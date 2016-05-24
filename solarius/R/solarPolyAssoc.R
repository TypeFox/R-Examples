# Run association analysis via polygenic model.
#
# @export
solarPolyAssoc <- function(formula, data, dir, kinship, 
  traits, covlist = "1",
  # input data to association 
  snpcovdata, snplist, snpind,
  assoc.options = "", assoc.settings = "", 
  assoc.df = 1, assoc.fix = FALSE, assoc.se = TRUE,
  # misc
  cores = getOption("cores"),
  ...,
  verbose = 0) 
{
  stopifnot(requireNamespace("tools"))

  tsolarPolyAssoc <- list()
  tsolarPolyAssoc$args <- proc.time()

  ### Step 1: Arguments
  mc <- match.call()
    
  if(!missing(snpcovdata)) {
    stopifnot(class(snpcovdata) == "matrix")
  }
  
  # check `snplist` / `snpind` format
  stopifnot(any(missing(snplist), missing(snpind)))
  
  if(missing(snplist)) {
    snplist <-  character()
  }
  if(missing(snpind)) {
    snpind <- integer()
  }

  if(!assoc.se) {
    assoc.settings <- c(assoc.settings, "option StandErr 0")
  }

  # cores
  if(is.null(cores)) {  
    cores <- 1
  }
    
  is.tmpdir <- missing(dir)
  
  ### step 2: process `snpcovdata`
  if(length(snplist)) {
    stopifnot(all(snplist %in% colnames(snpcovdata)))
    snpcovdata <- snpcovdata[, snplist, drop = FALSE]
  }
  if(length(snpind)) {
    stopifnot(all(snpind <= ncol(snpcovdata)))
    snpcovdata <- snpcovdata[, snpind, drop = FALSE]
  }
  
  snps <- colnames(snpcovdata)
  num.snps <- length(snps)
  snpcovdata <- data.frame(ID = rownames(snpcovdata), snpcovdata, stringsAsFactors = FALSE)
  
  fields <- matchIdNames(names(data))
  id.field <- names(fields)[1]

  names(snpcovdata)[1] <- id.field
  
  stopifnot(!any(names(snpcovdata)[-1] %in% names(data)))
  gdata <- merge(data, snpcovdata, by = id.field, all = TRUE)
  
  ### step 3: SOLAR dir
  if(is.tmpdir) {
    dir <- tempfile(pattern = "solarPolyAssoc-")
  }
  if(verbose > 1) cat(" * solarPolyAssoc: parameter `dir` is missing.\n")
  if(verbose > 1) cat("  -- temporary directory `", dir, "` used\n")

  ### Step 4: compute a polygenic model by calling `solarPolygenic`
  tsolarPolyAssoc$polygenic <- proc.time()

  if(!missing(formula)) {
    formula <- update(formula, paste0("~ . + ", snps[1], "()"))
  } else {
    covlist <- c(covlist, paste0(snps[1], "()"))
  }
  
  out <- solarPolygenic(formula, gdata, dir,
    kinship, traits, covlist, ..., verbose = verbose)

  out$assoc <- list(call = mc,
    cores = cores,
    snps = snps, num.snps = num.snps,
    # input par
    snplist = snplist, snpind = snpind,
    assoc.informat = "snpcovdata", assoc.outformat = "df",
    # SOLAR options/settings
    assoc.options = assoc.options, assoc.settings = assoc.settings, assoc.fix = assoc.fix, assoc.df = assoc.df)

  ### Step 4.1: add assoc.-specific slots to `out`
  tsolarPolyAssoc$preassoc <- proc.time()

  ### step 5: run assoc  
  tsolarPolyAssoc$runassoc <- proc.time()

  parallel <- (cores > 1)
  if(parallel) {
    # load required R package doParallel
    stopifnot(requireNamespace("doParallel", quietly = TRUE))
    
    doParallel::registerDoParallel(cores = cores)
  }

  if(!parallel) {
    snpf <- run_poly_assoc(out, dir, out$assoc$snps)
  } else {
    num.snps <- out$assoc$num.snps
    
    num.gr <- 10 * cores 
    gr <- cut(1:num.snps, breaks = seq(1, num.snps, length.out = num.gr + 1), include.lowest = TRUE)

    snpf <- llply(1:nlevels(gr), function(i) {
      # snps of a group `i`
      gr.i <- levels(gr)[i]
      snps.i <- snps[gr %in% gr.i]
      
      if(length(snps.i) == 0) {
        return(NULL)
      }

      if(out$verbose) cat(" * solarPolyAssoc: ", i, "/", nlevels(gr), "batch...\n")

      # copy `dir`      
      files.dir <- list.files(dir, include.dirs = TRUE, full.names = TRUE)
      dir.assoc <- tempfile(pattern = paste0("solarPolyAssoc-", i, "-"))
      
      stopifnot(dir.create(dir.assoc, showWarnings = FALSE, recursive = TRUE))
      stopifnot(file.copy(from = files.dir, to = dir.assoc, recursive = TRUE))
      
      # run
      snpf <- run_poly_assoc(out, dir.assoc, snps.i)
      
      # clean
      unlink(dir.assoc, recursive = TRUE)
      
      return(snpf)
    }, .parallel = parallel)

    #snpf <- llply(1:num.snps, function(i) {
    #  if(out$verbose) cat(" * solarPolyAssoc: ", i, "/", num.snps, "snp...\n")
    #
    #  # copy `dir`      
    #  files.dir <- list.files(dir, include.dirs = TRUE, full.names = TRUE)
    #  dir.assoc <- tempfile(pattern = paste0("solarPolyAssoc-", i, "-"))
    #  
    #  stopifnot(dir.create(dir.assoc, showWarnings = FALSE, recursive = TRUE))
    #  stopifnot(file.copy(from = files.dir, to = dir.assoc, recursive = TRUE))
    # 
    #  # run
    #  snpf <- run_poly_assoc(out, dir.assoc, snps[i])
    #  
    #  # clean
    #  unlink(dir.assoc, recursive = TRUE)
    #  
    #  return(snpf)
    #}, .parallel = parallel)
    
    ret <- try({
      rbindlist(snpf)
    })

    if(class(ret)[1] != "try-error") {
      snpf <- ret
    }
  }

  out$snpf <- snpf
  
#-------------
if(FALSE) {
  model.path <- paste0(paste(out$traits, collapse = "."), "/", tools::file_path_sans_ext(out$solar$model.filename))
  
  if(FALSE) {
    cmd.proc_def <- c(get_proc_poly_assoc(snps))
    cmd.proc_call <- c(paste0("poly_assoc ", "\"", model.path, "\""))
  
    cmd <- c(cmd.proc_def, cmd.proc_call)
    ret <- solar(cmd, dir, result = TRUE)
  } else {
    if(length(out$traits) == 1) {
      cmd.proc_def <- c(get_proc_poly_assoc2(snps, "sex*$snp", 2))
      cmd.proc_call <- c(paste0("poly_assoc2 ", "\"", model.path, "\""))
  
      cmd <- c(cmd.proc_def, cmd.proc_call)
      ret <- solar(cmd, dir, result = TRUE)
    } else if(length(out$traits) == 2) {
      cmd.proc_def <- c(get_proc_poly_assoc2(snps, paste0("${snp}_2(", out$traits[2], ")"), 2))
      cmd.proc_call <- c(paste0("poly_assoc2 ", "\"", model.path, "\""))
  
      cmd <- c(cmd.proc_def, cmd.proc_call)
      ret <- solar(cmd, dir, result = TRUE)
    }
  }
  
  snpf <- try({
    tab <- fread(file.path(dir, "out.poly.assoc"), sep = " ", header = FALSE)
    if(ncol(tab) == 4) {
      setnames(tab, c("SNP", "pSNP", "pSNPc", "pSNPi"))
    }
  })
  
  out$snpf <- snpf
} 

#-------------

  if(is.tmpdir) {
    unlink(dir, recursive = TRUE)
    if(verbose > 1) cat("  -- solarPolyAssoc: temporary directory `", dir, "` unlinked\n")
  }

  ### return
  tsolarPolyAssoc$return <- proc.time()
  out$assoc$tprofile$tproc$tsolarPolyAssoc <- tsolarPolyAssoc

  out$assoc$tprofile <- try({
    procc_tprofile(out$assoc$tprofile)
  }, silent = TRUE)
  
  oldClass(out) <- c("solarPolyAssoc", "solarAssoc", oldClass(out))
  return(out)
}

run_poly_assoc <- function(out, dir, snps)
{
  model.path <- paste0(paste(out$traits, collapse = "."), "/", tools::file_path_sans_ext(out$solar$model.filename))
  assoc.fix <- out$assoc$assoc.fix
  assoc.df <- out$assoc$assoc.df
  assoc.settings <- out$assoc$assoc.settings
  
  cmd.proc_def <- c(get_proc_poly_assoc(snps, df = assoc.df, assoc.fix = assoc.fix, assoc.settings = assoc.settings))
  cmd.proc_call <- c(paste0("poly_assoc ", "\"", model.path, "\""))
  
  cmd <- c(cmd.proc_def, cmd.proc_call)
  ret <- solar(cmd, dir, result = TRUE)

  snpf <- try({
    tab <- fread(file.path(dir, "out.poly.assoc"), sep = " ", header = FALSE)
    if(ncol(tab) == 3) {
      setnames(tab, c("SNP", "chi", "pSNP"))
    }
  })
}

get_proc_poly_assoc <- function(snps, df = 1, assoc.fix = FALSE, assoc.settings = "")
{
paste0("proc poly_assoc {model} {\
\
  # write to file\
	set f [open \"out.poly.assoc\" \"w\"]\
	set f_loglik [open \"out.poly.assoc.loglik\" \"w\"]\
	\
	foreach snp [list ", paste(snps, collapse = " "), "] {\
    # model 0\
	  load model $model\
	  set loglike_0 [loglike]\
    \
	  # model 1\
	  load model $model\
	  set pval_full \"NA\"\
	  \
	  covariate $snp\
	  ",
	  ifelse(assoc.fix,
	    "constraint delete_all\
	    \
	    set traits [trait]\
	    for {set i 0} {$i < 2} {incr i} {\
        set ti [lindex $traits $i]\
        fix \"h2r(${ti})\"\
        fix \"e2(${ti})\"\
        fix \"sd(${ti})\"\
      }\
      fix rhoe\
      fix rhog\
	    ",
	    "\
	    "),
	  assoc.settings,
	  "\
	  if {![catch {maximize -q}]} {\
      #save model snp\
  	  set loglike_1 [loglike]\
      \
	    set D [expr 2 * ($loglike_1 - $loglike_0)]\
	    if {$D < 0} {\
	      set pval_full 1\
	    } else {\
  	    set pval_full [chi -number $D ", df, "]\
  	  }\
	  }\
	  \
	  # put\
	  puts $f_loglik \"$snp $loglike_0 $loglike_1\"\
	  puts $f \"$snp $D $pval_full\"\
	}\
	\
	close $f\
	close $f_loglik\
}\
")
}

get_proc_poly_assoc2 <- function(snps, cov2, df2)
{
paste0("proc poly_assoc2 {model} {\
\
  # write to file\
	set f [open \"out.poly.assoc\" \"w\"]\
	set f_loglik [open \"out.poly.assoc.loglik\" \"w\"]\
	\
	foreach snp [list ", paste(snps, collapse = " "), "] {\
	  # extra variable\
	  define ${snp}_2 = $snp\
	  \
    # model 0\
	  load model $model\
	  set loglike_0 [loglike]\
    \
	  # model 1\
	  load model $model\
	  set pval_common \"NA\"\
	  \
	  covariate $snp\
	  if {![catch {maximize -q}]} {\
  	  set loglike_1 [loglike]\
	    \
	    set D [expr 2 * ($loglike_1 - $loglike_0)]\
	    if {$D < 0} {\
	      set pval_common 1\
	    } else {\
  	    set pval_common [chi -number $D 1]\
  	  }\
	  }\
	  \
	  # model 2\
	  load model $model\
	  set pval_interaction \"NA\"\
	  set pval_full \"NA\"\
	  \
	  covariate $snp ", cov2, "\
	  if {![catch {maximize -q}]} {\
  	  set loglike_2 [loglike]\
	    \
	    set D [expr 2 * ($loglike_2 - $loglike_1)]\
	    if {$D < 0} {\
	      set pval_interaction 1\
	    } else {\
  	    set pval_interaction [chi -number $D 1]\
  	  }\
      \
	    set D [expr 2 * ($loglike_2 - $loglike_0)]\
	    if {$D < 0} {\
	      set pval_full 1\
	    } else {\
  	    set pval_full [chi -number $D ", df2, "]\
  	  }\
	  }\
	  \
	  # put\
	  puts $f_loglik \"$snp $loglike_0 $loglike_1 $loglike_2\"\
	  puts $f \"$snp $pval_full $pval_common $pval_interaction\"\
	}\
	\
	close $f\
	close $f_loglik\
}\
")
}

#--------------------
# Class solarPolyAssoc
#--------------------

#' S3 class solarPolyAssoc.
#'
#' @name solarPolyAssoc
#' @rdname solarPolyAssoc
#'
#' @param x 
#'    An object of class \code{solarPolyAssoc}.
#' @param ...
#'    Additional arguments.
#'
#' @exportClass solarPolyAssoc


#' @rdname solarPolyAssocClass
#' @export
print.solarPolyAssoc <- function(x, ...)
{
  cat("\nCall: ")
  print(x$assoc$call)
  
  cat("\n Input SNP data:\n")
  switch(x$assoc$assoc.informat,
    "snpcovdata" = cat("  *  ", x$assoc$num.snps, " SNP covariates passed by `snpcovdata` argument\n", sep = ""),
    stop("switch error")
  )
  
  cat("\n Output results of association:\n")
  if(x$assoc$assoc.outformat == "df") {
    cat("\n  *  Table of association results (first 5 out of ", nrow(x$snpf), " rows):\n", sep = "")
    print(head(x$snpf, 5))
  }
  
  t1 <- with(x$assoc$tprofile$tprocf$tsolarPolyAssoc, elapsed.diff[unit == "runassoc"])
  t2 <- x$assoc$tprofile$cputime.sec
  cat("\n CPU time on ", x$assoc$cores, " core(s):\n", 
    " -- ", format(.POSIXct(t1, tz = "GMT"), "%H:%M:%S"), " -- association\n", 
    " -- ", format(.POSIXct(t2, tz = "GMT"), "%H:%M:%S"), " -- in total\n", 
    sep = "")
}
