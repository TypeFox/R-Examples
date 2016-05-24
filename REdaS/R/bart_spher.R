bart_spher <- function(x, use = c("everything", "all.obs", "complete.obs", "na.or.complete", "pairwise.complete.obs")){
  
  this_call <- match.call()
  use <- match.arg(use)
  warn <- FALSE
  
  if(any(is.na(x)) && (use %in% c("everything", "all.obs"))){
    stop(paste0('"x" contains missing values and use = "', use, '"\n  Check your data and set the "use" argument appropriately.'))
  }

  cormat <- tryCatch(cor(x = x, use = use), error = function(e) "error")
  
  if((cormat == "error") && (use == "complete.obs")) stop('"x" contains missing values and use = "complete.obs" produced an error. Check your data.')
  
  if(any(is.na(cormat))){
    if(use == "na.or.complete") stop()
    if(use == "pairwise.complete.obs") stop()
  }
  
  if(!any(is.na(x))){
    n <- nrow(x)
  } else if(use == "pairwise.complete.obs"){
    warn <- TRUE
    n_nNA <- which(rowSums(!is.na(x)) >= 2) # observations that are not NA (>= 2 responses)
    if(length(n_nNA) == 0) stop('no observations left to perform the test.')
    n <- sum(1.0 - is.na(x[n_nNA,])) / ncol(x[n_nNA,])
  } else {
    warn <- TRUE
    n <- nrow(na.omit(x))
  }
  k <- ncol(x)
  
  X2 <- - (n - 1 - (2 * k + 5) / 6) * determinant(cormat, logarithm = TRUE)[["modulus"]]
  df <- k * (k - 1) / 2
  p  <- pchisq(X2, df, lower.tail = FALSE)
  
  res <- structure(
    list("call"    = this_call,
         "x"       = x,
         "cormat"  = cormat,
         "use"     = use,
         "n"       = n,
         "k"       = k,
         "X2"      = as.numeric(X2), # remove attributes
         "df"      = df,
         "p.value" = as.numeric(p),
         "warn"    = warn),
    "class" = "bart_spher"
  )
  
  return(res)
  
}
