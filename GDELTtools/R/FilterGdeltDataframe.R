FilterGdeltDataframe <- function(x, 
                                 filter, 
                                 allow.wildcards=FALSE, 
                                 use.regex=FALSE,
                                 verbose=TRUE) {
  # 'or' within values for a field, 'and' across fields
  #
  # ex: FilterGdeltDataframe(my.df, list(ActionGeo_ADM1Code=c("NI", "US"), ActionGeo_CountryCode="US"))
  # This keeps rows with (ActionGeo_ADM1Code=NI AND ActionGeo_CountryCode=US) OR
  #   (ActionGeo_ADM1Code=US AND ActionGeo_CountryCode=US)
  if(verbose) cat("Filtering\n\n")
  if(use.regex) {
    filter.results <- laply(1:length(filter), function(fi) {
      field.results <- laply(.data=filter[[fi]], .fun=function(v) {
        grepl(v, x[names(filter)[fi]][,1])
      }, .drop=FALSE)
      if(is.array(field.results)) return(apply(field.results, 2, any))
      else return(field.results)
    })
  } else if(allow.wildcards) {
    filter.results <- laply(1:length(filter), function(fi) {
      field.results <- laply(.data=filter[[fi]], .fun=function(v) {
        v <- gsub("*", "[:alnum:]*", v, fixed=TRUE)
        grepl(v, x[names(filter)[fi]][,1])
      }, .drop=FALSE)
      if(is.array(field.results)) return(apply(field.results, 2, any))
      else return(field.results)
    })
  } else {
    filter.results <- laply(1:length(filter), function(fi) {
      field.results <- laply(.data=filter[[fi]], .fun=function(v) x[names(filter)[fi]]==v, .drop=FALSE)
      if(is.array(field.results)) return(apply(field.results, 2, any))
      else return(field.results)
    })
  }
  
  if(is.array(filter.results)) rows.to.keep <- apply(filter.results, 2, all)
  else rows.to.keep <- filter.results
  
  out <- x[rows.to.keep,]
  
  # remove NA values for filtered fields
  for(i in 1:length(filter)) {
    keep.rows <- !is.na(out[,names(filter)[i]])
    out <- out[keep.rows,]
  }
  
  return(out)
}
