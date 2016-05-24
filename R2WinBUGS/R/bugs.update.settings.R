"bugs.update.settings" <-
  function (n.burnin, bugs.directory)
{
  if(is.R()) {
    .fileCopy <- file.copy
    .regexpr <- function(...) regexpr(..., useBytes = TRUE)
    .sub <- function(...) sub(..., useBytes = TRUE)
    .writeBin <- 
        if(getRversion() >= "2.10") 
            function(...) writeBin(..., useBytes = TRUE)
        else 
            writeBin
  } else {
    .fileCopy <- splus.file.copy
    .regexpr <- regexpr
    .sub <- sub
    .writeBin <- writeBin
  }
  
  
  char.burnin <- as.character(n.burnin - 1)
  registry <- readBin(file.path(bugs.directory, "System/Rsrc/Registry.odc"),
                      "character", 400, size=1, endian="little")
  locale <- Sys.getlocale("LC_CTYPE")
  Sys.setlocale("LC_CTYPE", "C")
  .fileCopy(file.path(bugs.directory, "System/Rsrc/Registry.odc"),
            file.path(bugs.directory, "System/Rsrc/Registry_Rsave.odc"),
            overwrite = TRUE)
  info <- registry[.regexpr("Int", registry, fixed = TRUE) > 0]
  while(.regexpr("\r", info) > 0) {
    newline <- .regexpr("\r", info)
    info <- substring(info, newline + 1)
    line <- substring(info, 1, .regexpr("\r", info) - 1)
    if(.regexpr("AdaptivePhase", line) > 0) {
      numpos <- .regexpr("Int", line, fixed=TRUE) + 4
      num <- substring(line, numpos)
      if(as.numeric(num) > n.burnin) {
        blanks <- rep(" ", nchar(num, type="chars") - nchar(char.burnin, type="chars"))
        num.new <- paste(paste(blanks, collapse=""), char.burnin, sep="")
        line.new <- .sub(num, num.new, line)
        registry <- .sub(line, line.new, registry)
      }
    }
  }
  Sys.setlocale("LC_CTYPE", locale)
  .writeBin(registry,
           file.path(bugs.directory, "System/Rsrc/Registry.odc"), endian="little")
}


splus.file.copy <- function(from, to, overwrite=FALSE)
{
    if(!file.exists(from))
        stop("File: ", from, " does not exist")
    if(!overwrite && file.exists(to))
        stop("File: ", to, " already exists and overwrite is FALSE")
    n <- file.info(from)$size
    z <- writeBin(readBin(from, what="integer", size=1, n=n), to, size=1)
    invisible(z)
}
