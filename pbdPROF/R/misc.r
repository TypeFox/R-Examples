# convert NULLS in a list to 0
null_2_zero <- function(x)
{
  lapply(x, function(z) if(is.null(z)) 0 else z)
}

# convert parsed fpmpi list into a dataframe
parsed_fpmpi_2_df <- function(x)
{
  ret <- lapply(x, null_2_zero)
  ret <- t(matrix(unlist(ret), ncol=length(ret)))
  
  ret <- as.data.frame(ret, stringsAsFactors=F)
  
  nm <- names(x[[1L]])
  colnames(ret) <- nm
  
  # cast numeric columns appropriately
  numeric <- which(nm != "Routine")
  
  for (i in numeric){
    ret[, i] <- as.numeric(ret[, i])
  }
  
  return( ret )
}
