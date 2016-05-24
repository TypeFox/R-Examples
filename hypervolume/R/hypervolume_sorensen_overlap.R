hypervolume_sorensen_overlap <- function(hvlist)
{
  message("Note that hypervolume_sorensen_overlap definition has changed since version 1.2.")
  
  value <- 2 * hvlist@HVList$Intersection@Volume / (hvlist@HVList$HV1@Volume + hvlist@HVList$HV2@Volume)
  return(value)
}