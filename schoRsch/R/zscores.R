# -----------------------------------------------
# Compute z scores by condition. Input is 
# (a) a data frame witch specified columns
# for DV and Conditions or (b) A vector
# as DV and a list of conditions.
# -----------------------------------------------
zscores <- function(data,
                    factors=NaN,
                    dv = NaN) {

# Check input:
#
# Data frames
if ((is.data.frame(data)) &
    (!is.nan(dv))) {

  # Get DV column
  if (is.character(dv)) {
    which_column <- which(names(data) == dv)
  } else if (is.numeric(dv)) {
    which_column <- dv 
  }
  # Check factor definition
  if (length(factors)<=1) {
    if (is.nan(factors))    {
      factors = 0*c(1:length(data[[which_column]]))
    } else {
      a <- factors
      factors <- list(data[,a])[[1]]
    }
  } else if (TRUE) {
      a <- factors
      factors <- list(as.list(data[,a]))[[1]]
  }
  # Compute z-scores
  zout <- unsplit(lapply(split(data[[which_column]], factors), scale), factors)
  
# Vectors
} else if (is.list(factors)) {

  # Compute z-scores
  zout <- unsplit(lapply(split(data, factors), scale), factors)
  
} else if (is.nan(factors)) {

  # Compute z-scores
  zout <- scale(data)
  zout <- zout[,]
  
} else {
  
  print("Error - possibly wrong input specified?")  
  stop(call. = TRUE)
}

if (max(is.na(zout)) > 0) {
  print("Output contains NAs, possibly due to cells with only one data point.")
}

# Result
zout

}