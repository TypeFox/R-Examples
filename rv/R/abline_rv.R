# ========================================================================
# abline.rv  -  a + bx
# ========================================================================
# change this to
#   mapply.rv("abline", ..., MoreArgs=NULL)
#

abline.rv <- function (a = NULL, b = NULL, h = NULL, v = NULL, ...) { ## NEW
  if (! anyisrv(a, b, h, v)) {
    return(graphics:::abline(a=a, b=b, h=h, v=v, ...))
  }
  line.sample <- rvpar("line.sample")
  if (!is.numeric(line.sample)) {
    stop("rvpar('line.sample') is not a number")
  }
  args <- list(FUN=graphics:::abline, a=a, b=b, h=h, v=v)
  nulls <- sapply(args, is.null)
  nullArgs <- names(nulls)[nulls]
  MoreArgs <- list(...)
  args[nullArgs] <- NULL
  MoreArgs[nullArgs] <- list(NULL)
  args$SAMPLESIZE <- line.sample
  args$MoreArgs <- MoreArgs
  do.call(rvmapply, args=args)
  invisible()
}

