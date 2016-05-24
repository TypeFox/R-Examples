# Preparation for row-wise imputation
# Version:                           0.2
# Date:                       2015-05-17
# Author: F.M., some contributions: T.S.

rowimpPrep <- function(data,
                       ID=NULL,
                       verbose=TRUE)
{
	
  data          <- as.data.frame(data)
  nam_dat       <- names(data)
  key           <- NULL
  IDc           <- ID
  origindMatrix <- is.na(data)
  
  if (length(ID) > 2) stop("argument 'ID' contains more than two elements!")
  if (!is.null(ID)) {
  	if (is.character(IDc)) ID <- which(is.element(nam_dat,IDc))
  	key <- data.frame(data[ ,ID])
  	names(key) <- nam_dat[ID]
  	data <- data[ ,-ID] } else {key <- NULL}
  block <- list()

  ### Test for zero variance
  pos_novar <- apply((asNumericMatrix(data)),2,var,na.rm=T)==0
  ignore <- ignored_data <- NULL
  if(sum(pos_novar) > 0){
  	warning(paste("Variable",
  								which(pos_novar),
  								"contains no variation and will be ignored but added again after imputation. ")
  					)
  	ignore <- which(pos_novar)
  	ignored_data <- data[,ignore]
  	data <- data[,-ignore]
  }

  ##org.names <- var.names <- names(data) # variable names
  l    <- ncol(data) # no. variables
  miss <- function(x) {any(is.na(x)) }
  mvar <- apply(data, 2, miss)  #which columns have missings values?

  ## Stop, when no missing values were found
  if(sum(mvar)==0) stop(paste("No missing values found."))

  p.impvar <- p.impvar2 <- (1:l)[mvar] #position of columns with missing

  ## values
  ## constant regressors?

  obscol <- l - length(p.impvar)
  
  # completely observed variables
  if(obscol == 1){
    comp <-  data[-p.impvar] 
  } else {
    comp <- data[ ,-p.impvar] 
  }


  ## values
  comp.names <- names(comp) # names of completely observed variables

  dat.isna   <- is.na(data) # indicator matrix for missing values
  co1 <- 0
  for (i in 1:length(p.impvar))
    ## creates list with positions for identical columns (blocks)
    {
      co1 <- co1+1
      if(length(p.impvar2) > 1)
        {
          aa <- (dat.isna[,p.impvar2[1]] == dat.isna[,p.impvar2]) - 1
          identical <- which(apply(aa, 2, sum) == 0)
          ## which columns are identical to the first column with missing
          ## values
          block[[co1]] <- p.impvar2[identical]
        } else if (length(p.impvar2) == 1){
          block[[co1]] <- p.impvar2
          break  # break if one column is left
        } else {
          break  # break if no column is left
        }
      p.impvar2 <- p.impvar2[-identical]
    }
  block.names <- lapply(block,FUN=function (x) names(data[x]))
  names(block.names) <- names(block) <- paste("block",1:length(block))
  if (verbose) {
    cat("variables with missing values ordered by identical patterns:\n")
    print(block.names)
  }
  
  x <- list("data"          = data,
  					"key"           = key,
  					"blocks"        = block,
  					"blockNames"    = block.names,
  					"compNames"     = comp.names,
  					"ignore"        = ignore,
  					"ignored_data"  = ignored_data,
  					"indMatrix"     = origindMatrix)
  
  class(x) <- "impprep"
  return(x)
}
