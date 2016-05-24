AssignSameRustlook <-
function(table,rustlook){
  # Assigns the same outlook ('rustlook') to all unpublished studies in the data set.
  # The assigned rustlook will override any previously assigned outlook. 
  #
  # Args: 
  #   table: The data set.
  #   rustlook: The outlook to be assigned to all unpublished studies.
  #
  # Returns: The data set with the rustlook assigned to all unpublished studies.

  # Examples:
  #   AssignSameRustlook(Hpylori,"negative")
  #   AssignSameRustlook(Hpylori,"current effect")

  # Dependencies:
  #   Calls: CompleteOutlooksFactor()

  table <- CompleteOutlooksFactor(table)
  table$outlook[which(table$outlook != "published")] <- rustlook
  return(table)
}
