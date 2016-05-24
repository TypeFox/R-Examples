ReturnIconColour <-
function(iIcon)
{
  iColour <- iIcon %% 9
  
  if (iColour == 1)
    IconColour <- "cc00ff55"
  if (iColour == 2)
    IconColour <- "cc00ffbf"
  if (iColour == 3)
    IconColour <- "cce65c00"
  if (iColour == 4)
    IconColour <- "ccff57ff"
  if (iColour == 5)
    IconColour <- "cc0066e3"
  if (iColour == 6)
    IconColour <- "cc000080"
  if (iColour == 7)
    IconColour <- "ccffc500"
  if (iColour == 8)
    IconColour <- "cc00c8ff"
  if (iColour == 0)
    IconColour <- "ccff008c"
  
  IconColour
}
