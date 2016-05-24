CompleteOutlooksFactor <-
function(table){
  # When importing a data set, R will include in the 'outlooks' vector only those outlooks found in the data set. 
  # This function completes the 'outlooks' vector with other outlooks not found in the data set,
  # so that any outlook may be summoned when imputing unpublished studies with a different outlook. 
  #
  # Args: 
  #   table: The data set with an 'outlook' column. 
  #
  # Returns: The same data set with the factor 'outlook' complete.

  # Example:
  #   greentea$outlook
  #   > Levels: no effect positive published
  #   greentea <- CompleteOutlooksFactor(greentea)
  #   greentea$outlook
  #   > 11 Levels: no effect positive published very positive ... very negative CL 

  # Dependencies:
  #   Callers: AssignSameRustlook()

  outlooks <- c("very positive", "positive", 
                "no effect", 
                "negative", "very negative", 
                "very positive CL", "positive CL", 
                "current effect", 
                "negative CL", "very negative CL")
  levels(table$outlook) <- union(levels(table$outlook),outlooks)
  return(table)
}
