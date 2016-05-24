### wideByFactor.R
###------------------------------------------------------------------------
### What: Reshape by factor levels - code
### $Id$
### Time-stamp: <2008-12-30 22:17:32 ggorjan>
###------------------------------------------------------------------------

wideByFactor <- function(x, factor, common, sort=TRUE, keepFactor=TRUE)
{
  ## --- Setup ---
  if(!is.data.frame(x)) stop("'x' must be a data frame")
  if(length(factor) != 1) stop("'factor' can be only of length one") 
  if(!is.factor(x[[factor]])) stop("column defined in 'factor' must be a factor")
  if(sort) x <- x[order(x[[factor]]), ]

  ## --- Extend by factors levels ---
  y <- x[common]
  if(keepFactor) y[factor] <- x[factor]
  levs <- levels(x[[factor]])

  ## Remove common and factor from the list of column names
  other <- names(x)
  other <- other[!(other %in% common) & !(other %in% factor)]

  ## Add all other columns but as a set for each level of a factor
  for(level in levs) {
    for(col in other) {
      ## add a column col
      y[paste(col, level, sep=".")] <- x[col]
      ## fill with NA for other levels than level
      y[x[factor] != level, paste(col, level, sep=".")] <- NA
      ## This filling migth be inefficient if there is large number
      ## of levels, since there will be quite a lot of filling.
    }
  }
  y
}

###------------------------------------------------------------------------
### wideByFactor.R ends here