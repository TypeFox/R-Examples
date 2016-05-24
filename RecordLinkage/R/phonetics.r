# phonetic.r: functions for phonetic coding

pho_h <- function(str)
{
  # check type
  if (typeof(str) != "character" && class(str) != "factor")
     stop(sprintf("Illegal data type: %s", typeof(str)))
  if (class(str) == "factor")
    str=as.character(str)
   out <- .C("pho_h", str, ans=character(length(str)),length(str),
             PACKAGE="RecordLinkage")
   if (any(is.na(str)))
    out$ans[is.na(str)]=NA
   dim(out$ans)=dim(str)
   dimnames(out$ans)=dimnames(str)
   return(out$ans)
}

soundex <- function(str)
{
  # check type
  if (typeof(str) != "character" && class(str) != "factor")
     stop(sprintf("Illegal data type: %s", typeof(str)))
  if (class(str) == "factor")
    str=as.character(str)
   out <- .C("soundex", as.character(str), ans=character(length(str)),length(str),
             PACKAGE="RecordLinkage")
   if (any(is.na(str)))
    out$ans[is.na(str)]=NA
   dim(out$ans)=dim(str)
   dimnames(out$ans)=dimnames(str)
   return(out$ans)
}
