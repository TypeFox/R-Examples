GenerateTabs <-
function(iNumberOfTabs)
{
  # 
  # Chr(9) is the horizontal tab character
  # Use 8 spaces instead of a tab character.
  
  sResult <- ""
  for (i in 1:iNumberOfTabs)
    sResult <- paste(sResult,"        ",sep="")
  
  sResult
}
