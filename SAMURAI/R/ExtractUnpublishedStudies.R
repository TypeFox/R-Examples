ExtractUnpublishedStudies <-
function(table, colname="outlook"){
  # Extracts the unpublished studies (outlook != "published") from a prepared data set.
  #
  # Args: 
  #   table: The data set, in table form. 
  #   colname: The column containing the outlooks for the studies in the data set.
  #
  # Returns: The subset of unpublished studies. 

  out <- table[which(with(table, get(colname)) != "published" ), ]
  ## Rename column as 'outlook'
  names(out)[which(names(out) == colname)] <- "outlook" 
  return(out)
}
