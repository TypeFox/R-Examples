lrvar <- function(x, type = c("Andrews", "Newey-West"), prewhite = TRUE, adjust = TRUE, ...)
{
  ## data and model
  x <- coredata(x)
  m <- lm(x ~ 1)
 
  ## vcov matrix
  type <- match.arg(gsub("-", "", tolower(type), fixed = TRUE),
    c("andrews", "neweywest"))
  vc <- switch(type,
    "andrews" = kernHAC(m, prewhite = prewhite, adjust = adjust, ...),
    "neweywest" = NeweyWest(m, prewhite = prewhite, adjust = adjust, ...))

  ## set names or drop dimension for univariate series
  if(NCOL(x) > 1) {
    rownames(vc) <- colnames(vc) <- colnames(x)
  } else {
    vc <- drop(vc)
  }
  return(vc)
}
