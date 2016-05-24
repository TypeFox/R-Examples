`dimnames.mefa` <-
function (x)
{
  out <- list(samp = rownames(x$xtab),
       taxa = colnames(x$xtab),
       segm = names(x$segm))
  if (is.null(out$segm))
       out$segm <- "undefined"
  return(out)
}

