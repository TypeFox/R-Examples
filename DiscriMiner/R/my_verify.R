my_verify <-
function(x, y, qualitative=FALSE, na.rm=na.rm)
{
  # x: matrix or data frame with explanatory variables
  # y: vector or factor with group memberships
  # qualitative: logical indicating verification for disqual
  # na.rm: logical indicating missing values in x
  
  # x matrix or data.frame
  if (is.null(dim(x))) 
    stop("\n'variables' is not a matrix or data.frame")
  # no missing values allowed when na.rm=FALSE
  if (!na.rm) {
    if (length(complete.cases(x)) != nrow(x))
      stop("\nSorry, no missing values allowed in 'variables'")    
  }
  # check lengths of x and y
  if (nrow(x) != length(y))
    stop("\n'variables' and 'group' have different lengths")
  # y vector or factor
  if (!is.vector(y) && !is.factor(y))
    stop("\n'group' must be a factor")
  # make sure y is a factor
  if (!is.factor(y)) y = as.factor(y)
  # no missing values in y
  if (any(!is.finite(y)))
    stop("\nNo missing values allowed in 'group'")
  # quantitative or qualitative variables?
  if (!qualitative)
  { # quantitative data
    # make sure is matrix
    if (!is.matrix(x)) x <- as.matrix(x)
    # only numeric values
    if (!is.numeric(x))
      stop("\n'variables' must contain only numeric values")    
  } else { # data frame with qualitative data
    # variables as data frame with factors
    fac_check = sapply(x, class)
    if (!is.data.frame(x) && any(fac_check != "factor"))
      stop("\nA data frame with factors was expected")
  }
  
  # verified inputs
  if (is.null(colnames(x))) 
    colnames(x) = paste(rep("X", ncol(x)), seq_len(ncol(x)), sep='') 
  if (is.null(rownames(x))) 
    rownames(x) = 1L:nrow(x)
  list(X=x, y=y)
}
