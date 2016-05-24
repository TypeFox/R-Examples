multicomp.reverse <- function(y, estimate.sign=1, ...) {
  mca <- y
  if (estimate.sign == 0) return(mca)   ## no change, return argument
  if (estimate.sign > 0)
    tml <- (mca$table[,"estimate"] < 0) ## locate negatives to be changed
  else
    tml <- (mca$table[,"estimate"] > 0) ## locate positives to be changed
  mca$table[tml,"estimate"] <- -mca$table[tml,"estimate",drop=FALSE]
  mca.table.tml..upper      <-  mca$table[tml,"upper",drop=FALSE]
  mca$table[tml,"upper"]    <- -mca$table[tml,"lower",drop=FALSE]
  mca$table[tml,"lower"]    <- -mca.table.tml..upper
  mca$lmat[,tml]            <- -mca$lmat[,tml,drop=FALSE]

  contrast.names <- dimnames(mca$table)[[1]][tml]
  minus.position <- regexpr("-",contrast.names)
  if (any(minus.position==-1))
    cat("mmc: At least one reversed contrast name did not have a '-' sign.\n     We appended a '-' sign.\n")
  ## message("mmc: At least one reversed contrast name did not have a '-' sign.\n     We appended a '-' sign.")
    ## warning("At least one reversed contrast name did not have a '-' sign.  We appended a '-' sign.")
  last.position  <- nchar(contrast.names)
  contrast.names.rev <-
    paste(substring(contrast.names, minus.position+1, last.position),
          substring(contrast.names, 1, minus.position-1),
          sep="-")

  dimnames(mca$table)[[1]][tml]    <- contrast.names.rev
  dimnames(mca$lmat)[[2]][tml]     <- contrast.names.rev
  if (!is.null(mca$height))
    names(mca$height)[tml] <- contrast.names.rev
  mca
}
