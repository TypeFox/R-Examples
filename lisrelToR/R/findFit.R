findFit <- function(Out,start)
{
  ### Extract fit statistics:
  IndStart <- start
  # Find end:
  IndEnd <- IndStart
  repeat
  {
    IndEnd <- IndEnd + 1
    if (!(grepl("\\s*",Out[IndEnd]) | grepl("=",Out[IndEnd]))) break
  }
  modTxt <- Out[IndStart:IndEnd]
  modTxt <- modTxt[grepl("=",modTxt)]
  modTxt <- strsplit(modTxt,split="=")
  modTxt <- lapply(modTxt,gsub,pattern="^\\s*",replacement="")
  modTxt <- lapply(modTxt,gsub,pattern="\\s*$",replacement="")
  motTxt <- lapply(modTxt,function(x)c(x[1],paste(x[-1],collapse="=")))
  
  return(data.frame(Statstic=sapply(modTxt,"[",1),Value=sapply(modTxt,"[",2)))
}