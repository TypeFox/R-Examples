
## experimental / probably not very efficient 
## function for converting a column of multinominal values into logical matrix
## x: SPC
## v: site-level variable name, must be a factor
multinominal2logical <- function(x, v) {
  
  if(class(x) != 'SoilProfileCollection')
    stop('`x` must be a SoilProfileCollection', call. = FALSE)
  
  # get internal ID name
  id <- idname(x)
  
  # site data only
  s <- site(x)
  
  if(class(s[[v]]) != 'factor')
    stop('`v` must be a factor', call. = FALSE)
  
  # construct formula for dcast
  fm <- paste0(idname(x), ' ~ ', v)
  
  # cast to wide format, filling non-NA entries with ID
  w <- dcast(s, fm, value.var=id, drop=FALSE) 
  
  # not done yet, neet to convert into logical
  # first column is the ID
  w.logical <- sapply(w[, -1], function(i) ! is.na(i))
  
  # merge ID back in
  w.final <- data.frame(w[, 1], w.logical[, levels(s[[v]])], stringsAsFactors = FALSE)
  names(w.final)[1] <- id
  
  return(w.final)
}

