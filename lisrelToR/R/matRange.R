matRange <- function(start,lisName,Out)
{
  if (grepl("IN THE FOLLOWING GROUP",Out[start])) return(c(start,start))
  
  if (missing(lisName)) lisName <- gsub("(^\\s*)|(\\s*$)","",Out[start])
  # Find end:
  IndEnd <- start
  repeat
  {
    IndEnd <- IndEnd + 1
    if (length(Out)==IndEnd) break
    if (Out[IndEnd]=="" & !grepl(lisName,Out[IndEnd+1]) & !grepl(lisName,Out[IndEnd-1]) & !grepl(" Note:",Out[IndEnd-1],ignore.case=TRUE)) break
  }
  return(c(start,IndEnd))
}