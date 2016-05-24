ExtractPublishedStudies <-
function(table, colname="outlook"){  
  # Extracts the published studies (outlook == "published") from a prepared data set.
  #
  # Args: 
  #   table: The data set, in table form. 
  #   colname: The column containing the outlooks for the studies in the data set.
  #
  # Returns: The subset of published studies. 

  # Note:
  #   Use with(get()) to treat column names (strings) as variables:
  #   The with() statement looks for variables inside a dataframe.
  #   The get() statement accesses the object 'colname'.

  out <- table[which(with(table, get(colname)) == "published" ), ]
  # Rename column as 'outlook'
  names(out)[which(names(out) == colname)] <- "outlook" 
  return(out)  
}
