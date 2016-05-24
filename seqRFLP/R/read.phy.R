read.phy <-
function(fil = NULL)
{
  if(is.null(fil)){
    stop("You have to specify the input phylip file.")
  }
  fil <- readLines(fil)
  fil = fil[grepl("[A-Za-z1-9]", fil)]
  if(length(fil) <= 1){
    stop("The input file contains only one row, is it in phy format?")
  }
  class(fil) <- "phy"
  return(fil)
}

