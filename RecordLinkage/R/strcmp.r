# strcmp.r: functions for string comparison


jarowinkler <- function(str1, str2, W_1=1/3, W_2=1/3, W_3=1/3, r=0.5)
{
  # check type
  if (typeof(str1) != "character" && class(str1) != "factor")
     stop(sprintf("Illegal data type: %s", typeof(str1)))
  if (class(str1) == "factor")
    str=as.character(str1)

  if (typeof(str2) != "character" && class(str2) != "factor")
     stop(sprintf("Illegal data type: %s", typeof(str2)))
  if (class(str2) == "factor")
    str=as.character(str2)

   if ((is.array(str1) || is.array(str2)) && dim(str1)!=dim(str2))
      stop ("non-conformable arrays")
   if(length(str1)==0 || length(str2)==0) return(double(0))
    ans <- .Call("jarowinklerCALL", str1, str2, #l1,l2,
             as.double(W_1), as.double(W_2), as.double(W_3), as.double(r),
             PACKAGE="RecordLinkage")

   if (any(is.na(str1),is.na(str2)))
      ans[is.na(str1) | is.na(str2)]=NA

   if (is.array(str1))
      return(array(ans,dim(str1)))

   if (is.array(str2))
      return(array(ans,dim(str2)))

   return(ans)
}

levenshteinDist <- function(str1, str2)
{
  # check type
  if (typeof(str1) != "character" && class(str1) != "factor")
     stop(sprintf("Illegal data type: %s", typeof(str1)))
  if (class(str1) == "factor")
    str=as.character(str1)

  if (typeof(str2) != "character" && class(str2) != "factor")
     stop(sprintf("Illegal data type: %s", typeof(str2)))
  if (class(str2) == "factor")
    str=as.character(str2)
   if ((is.array(str1) || is.array(str2)) && dim(str1)!=dim(str2))
      stop ("non-conformable arrays")
   if(length(str1)==0 || length(str2)==0) return(integer(0))
   l1 <- length(str1)
   l2 <- length(str2)
   out <- .C("levenshtein", as.character(str1), as.character(str2),l1,l2,
             ans=integer(max(l1,l2)), 
             PACKAGE="RecordLinkage")

   if (any(is.na(str1),is.na(str2)))
      out$ans[is.na(str1) | is.na(str2)]=NA

   if (is.array(str1))
      return(array(out$ans,dim(str1)))

   if (is.array(str2))
      return(array(out$ans,dim(str2)))
      
   return(out$ans)
}

levenshteinSim <- function(str1, str2)
{
	return (1-(levenshteinDist(str1,str2)/pmax(nchar(str1),nchar(str2))))
}
