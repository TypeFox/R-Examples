set.lv <- function (object, lv, verbose = TRUE)
{
  
  iter <- get.iter (object)
  N <- nrow (object$scores)
  totlv <- dim (object$coefficients) [1]
  maxnlv <- ceiling (N / 10)

  if (lv == 'crop')
  {
    if (lv > maxnlv) new.lv <- maxnlv
    else new.lv <- lv
  }
  else
  {
    if (lv > totlv) stop (paste ('LV available:', totlv))
    else new.lv <- lv
  }
  names (new.lv) <- NULL
  new.lv <- as.integer (new.lv)
   
  ## $metapls object
  part2 <- object$metapls
  part2$current.iter <- iter
  part2$current.lv <- new.lv
  object$metapls <- part2

  ## Report
  if (verbose) print (object)

  ## Output
  invisible (object) 
}

