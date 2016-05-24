FillInMissingEffectSizeDefaults <-
function(exhibit, vault){
  # Replace missing values in the 'exhibit' vector 
  # with corresponding entries in the 'vault' vector.
  #
  # Args: 
  #   exhibit: the vector of effect sizes to update.
  #   vault: the vector of effect sizes to draw from.
  #
  # Returns:  The exhibit vector with missing values replaced by 
  #           corresponding entries in the vault vector. 

  # Example: 
  #   u <- c(2,NA)
  #   c <- c(5,3)
  #   FillInMissingEffectSizeDefaults(c,u) 
  #   > 2 3

  # Note: 
  #  There is no easy way to have elementwise addition of vectors while ignoring NA's.
  #  Setting NA's to zeros helped workaround this issue.

  vault.miss <- is.na(vault) ## missingness indicators
  retain <- exhibit * vault.miss
  vault[which(is.na(vault))] <- 0 ## change NA to zeros
  out <- vault + retain
  return(out)
}
