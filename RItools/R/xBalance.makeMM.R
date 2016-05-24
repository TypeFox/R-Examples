xBalance.makeMM <- function(tfm, dat) {
  attr(tfm, "intercept") <- 1
  mf <- model.frame(tfm, dat, na.action=na.pass) ##na.pass leaves the NAs in dat in the model.frame
  tlbl <- names(mf)
  names(tlbl) <- as.character(tlbl)
  clist <- lapply(mf,
                  function(x) {
                    if (is.factor(x))
                      structure(diag(nlevels(x)),
                                dimnames=list(levels(x), levels(x)))
                    else NULL })
  clist <- clist[!sapply(clist, is.null)]
  mm <- model.matrix(tfm,mf,contrasts.arg=clist)

  # creat a little look up table of the original variables
  assign <- attr(mm, "assign")[-1] # these will match up to the original variables
  originals <- attr(terms(tfm), "term.labels")[assign]

  tmp <- mm[, -1, drop = FALSE]
  attr(tmp, "originals") <- originals
  return(tmp)
}
