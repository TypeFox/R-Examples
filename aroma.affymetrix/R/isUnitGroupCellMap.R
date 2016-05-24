setMethodS3("isUnitGroupCellMap", "matrix", function(this, ...) {
  fields <- c("unit", "group", "cell");

  # Got the columns?
  if (!all(fields %in% colnames(this)))
    return(FALSE);

  # Are they numeric?
  if (!is.numeric(this))
    return(FALSE);

  TRUE;
})


setMethodS3("isUnitGroupCellMap", "data.frame", function(this, ...) {
  fields <- c("unit", "group", "cell");

  # Got the columns?
  if (!all(fields %in% colnames(this)))
    return(FALSE);

  # Are they numeric?
  for (field in fields) {
    if (!is.numeric(this[[field]]))
      return(FALSE);
  }
  
  TRUE;
})


setMethodS3("isUnitGroupCellMap", "default", function(this, ...) {
  FALSE;
})


############################################################################
# HISTORY:
# 2008-03-11
# o Created.
############################################################################
