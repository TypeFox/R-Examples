read.nxs <-
function(fil = NULL)
{
  if(is.null(fil))
  {stop("You have to specify the input nexus file.")}
  fil2 <- readLines(fil)
  fil3 = fil2[grepl("[A-Za-z0-9\\[]", fil2)]
  if(any(grepl("^[\\[]", fil3)))
  {stop("Currently no anotations are not permitted.\n Please delete the contents between \\[ and \\]\n.")}
  class(fil3) <- "nxs"
  return(fil3)
}

