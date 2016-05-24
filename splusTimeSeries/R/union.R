"unionPositions" <- 
function(..., localzone = FALSE, matchtol = 0)
{
  ## Make a union of numeric or calendar position objects using
  ## localzone and matchtol as in the series merge and align funs
  if(matchtol < 0) stop("matchtol must be >= 0")
  posunion2 <- function(pos1, pos2, localzone, matchtol)
    {
      ## find the union of two positions objects with zone and tolerance
      len1 <- length(pos1)
      len2 <- length(pos2)
      if(!len1)
        return(pos2)
      if(!len2)
        return(pos1)
      if(identical(pos1, pos2))
        return(pos1)
      if(localzone) {
        savezone <- pos1@time.zone
        savezone2 <- pos2@time.zone
        saveformat <- pos1@format
        ## only easy way to convert so comparison is done in local zone is
        ## to go to character (w/out time zone) and read back as default zone
        ## since timeZoneConvert goes the other direction
        pos1@format <- "%m %d %Y %H %M %S %N"
        pos2@format <- "%m %d %Y %H %M %S %N"
        pos1 <- timeDate(as(pos1, "character"), in.format = 
                         "%m %d %Y %H %M %S %N")
        pos2 <- timeDate(as(pos2, "character"), in.format = 
                         "%m %d %Y %H %M %S %N")
      }
      ## put all positions together, then remove sufficiently close
      ## duplicates
      pos1 <- sort(c(pos1, pos2))
      ok <- abs(diff(pos1)) > matchtol
      if(length(ok) && !all(ok)) {
        ## make it the right length
        ok <- c(ok, TRUE)
        pos1 <- pos1[ok]
      }
      if(localzone) {
        pos1@format <- saveformat
        pos1 <- timeZoneConvert(pos1, savezone)
      }
      pos1
    }
  arglist <- list(...)
  if(!length(arglist))
    return(numeric(0))
  pos <- sort(arglist[[1]])
  arglist <- arglist[-1]
  for(arg in arglist)
    pos <- posunion2(pos, sort(arg), localzone, matchtol)
  pos
}

