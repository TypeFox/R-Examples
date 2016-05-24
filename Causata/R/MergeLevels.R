MergeLevels <- function(this, ...) {
  # generic function to merge levels in factors
  UseMethod("MergeLevels", this)
}

# Takes a factor (nominal) variable and merges the least populated levels into 'Other' until
# we are left with only maxlevels levels
###############################################################################

MergeLevels.factor <- function(this, max.levels, other.name="Other", ...) {
  
  #Check that the arguments are valid
  if (max.levels < 2){
    stop("Error: cannot merge to fewer than 2 levels")
  }
  if (!is.factor(this)){
    stop("Column is not a factor - cannot merge")
  }
  if (length(levels(this)) <= max.levels){
    #cat("Factor only has",length(levels(column)),"levels - nothing to do \n")
    return(this)
  }
  
  # Get list of levels to merge into other, this works by getting a sorted table of levels with
  # largest levels first.  Levels at the back of the list will be merged.
  levelsToMerge <- names(sort(table(this), decreasing=TRUE)[ max.levels : length(levels(this)) ])
  
  # get the current levels before merge
  levelsOld <- levels(this)
  
  # get an index of levels to merge
  idx <- levelsOld %in% levelsToMerge
  
  # create a new vector of levels where smallest levels are merged
  levelsNew <- levelsOld
  levelsNew[idx] <- other.name
  
  # merge smallest levels
  levels(this) <- levelsNew
  
  return(this)
}
