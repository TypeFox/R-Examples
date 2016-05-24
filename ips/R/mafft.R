## PACKAGE: ips
## CALLED BY: USER, megaptera
## AUTHOR: Christoph Heibl (at gmx.net)
## LAST UPDATE: 2014-07-30

mafft <- function(x, y, add, method = "auto", maxiterate = 0, 
                  op = 1.53, ep = 0.0, gt, options, path, quiet){
  
  ## CHECKS and DEFINITIONS
  ## ----------------------
  if ( !inherits(x, "DNAbin") ) stop("'x' is not of class 'DNAbin'")
  os <- .Platform$OS
  if ( missing(quiet) ) quiet <- TRUE
  qut <- ifelse(quiet, " --quiet ", " ")
  if ( missing(path) ) path <- "/usr/local/bin/mafft"
  
  ## temporary input/output files
  ## ----------------------------
#   fns <- c("mafftin.fas", "mafftpro.fas", "mafftout.fas")
  fns <- vector(length = 3)
  for ( i in seq_along(fns) ) 
    fns[i] <- tempfile(pattern = "mafft", tmpdir = tempdir(), fileext = ".fas")
  unlink(fns[file.exists(fns)])
  method <- match.arg(method, c("auto", "localpair", "globalpair", 
                                "genafpair", "parttree",
                                "retree 1", "retree 2"))
  
  ## guide tree
  ## ----------
  if ( missing(gt) ){
    gt <- ""
  } else {
    if (!inherits(gt, "phylo"))                                          
      stop("object \"gt\" is not of class \"phylo\"")
    if ( !all(names(x) %in% gt$tip.label) )
      stop("guide tree does not match sequence names")
    gt$tip.label <- match(names(x), gt$tip.label)
    if ( !is.binary.tree(gt) )
      gt <- multi2di(gt)
    if ( is.null(gt$edge.length) )
      gt$edge.length <- rep(1, nrow(gt$edge))
    phylo2mafft(gt)
    gt <- " --treein tree.mafft "
  }
  
  ## additional arguments specified by the user
  ## ------------------------------------------
  if ( missing(options) ){
    options <- ""
  } else {
    options <- match.arg(options, c("--adjustdirection",
                                    "--adjustdirectionaccurately"))
    options <- paste(options, collapse = " ")
  }
  
  ## write input files and prepare call to MAFFT 
  ## -------------------------------------------
  if ( missing(y) ){
    write.fas(x, fns[1])
    call.mafft <- paste(path, " --", method, " --",   	
                        "maxiterate ", maxiterate, qut, "--op ", op, 		
                        " --ep ", ep, gt, " ", options, " ", 
                        fns[1], " > ", fns[3], sep = "")
  }
	else {
	  if ( !inherits(y, "DNAbin") ) stop("'y' is not of class 'DNAbin'")
    if ( missing(add) ) add <- "addprofile"
    add <- match.arg(add, c("add", "addprofile"))
    add <- paste("--", add, sep = "")
	  write.fas(x, fns[1])
	  write.fas(y, fns[2])
	  call.mafft <- paste(path, qut, add, fns[2], fns[1], ">", fns[3])
	}
  if ( !quiet ) message(call.mafft)
  
  ## execute MAFFT on UNIX
  ## ---------------------
  if ( os == "unix" ){
    system(call.mafft, intern = FALSE, ignore.stdout = FALSE)
    res <- length(scan(fns[3], what = "c", quiet = TRUE))
    if ( res != 0 ) res <- read.fas(fns[3])
  }
  ## execute MAFFT on WINDOWS
  ## ------------------------
  else {
    res <- system(call.mafft, intern = TRUE, ignore.stderr = FALSE)
    if ( length(grep("error|ERROR", res)) > 0 ){
      res <- 0
    }
    else {
      res <- read.fas(fns[3])
    }
  }
  unlink(fns[file.exists(fns)])
  return(res)
}

