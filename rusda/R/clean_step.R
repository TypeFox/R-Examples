##                       clean_step                  ##
##     This code is part of the rusda package        ##
##     F.-S. Krah 2015 (last update: 2015-07-11)     ##

clean_step <- function(x, syns, taxa, spec_type, synonyms_incl)
{
  if(is.null(x)) { x <- "nodata" }
  if(spec_type == "plant" & x[[1]] == "nodata" & length(x) == 1) { x }
  else{
  regex <- c(":\\s*-")
  del <- paste(regex, collapse = "|")
  if(length(grep(del, x)) > 0)   x <- x[-grep(del, x)]
  
  # delete empty strings
  if(length(which(nchar(x) == 0)) > 0) x <- x[-which(nchar(x) == 0)] 
  
  # flush left
  regex <- c(paste(intToUtf8(226),"\u0089", intToUtf8(161), sep=""),
             paste(intToUtf8(194),"\u0089",intToUtf8(161), sep=""), 
             paste(intToUtf8(194),"\\s*",sep=""),
             paste("=",intToUtf8(194),sep=""),"=", "Variant\\sspelling", 
             paste(intToUtf8(195),"\u0083\u0097",sep=""),":",
             paste(intToUtf8(195),"\u0083",sep=""), "\u0097",
             "BPI\\s[0-9]*", ",")
  del <- paste(regex, collapse = "|")
  if(length(grep(del, x)) > 0) x <- gsub(del, "", x)
  
  if(length(grep("^\\s+|\\s+$",x)) > 0) x <- gsub("^\\s+|\\s+$", "", x)
  if(length(grep("f. sp.", x)) > 0) x <- gsub("f. sp.", "f.sp.", x)
  
  # select first 4 words
  x <- lapply(x, function(x){
    con <- grep(paste(c("f\\.","var\\.","subsp\\.","f.\\ssp."), collapse="|" ), x)
    if(length(con) > 0){x <- word(x, 1, 4)}
    if(length(con) == 0 & sapply(gregexpr("[A-z]\\W+", x), length) >= 2){x <- word(x, 1, 2)}
    return(x)
  })
  
  # delete rows with taxa names or synonyms
  if(synonyms_incl==TRUE)
  {
    if(length(grep(taxa[[1]], x)) > 0) x <- x[-grep(taxa[[1]], x)]
    syns <- paste(syns[[1]], collapse="|")
    if(length(grep(syns[[1]], x)) > 0) x <- x[-grep(syns[[1]], x)]
  }

  # delete rows matching regex
  regex <- c("NA", "litter,\\sdecayed", "Substrate\\sUndetermined",
             "leaves,\\sdecayed", "water,\\sfresh","foam\\sHong", "root,\\sdecayed", 
             "bark,\\sdecayed", "twig\\sVenezuela", "humus\\sMexico","-\\s,",
             "Japan\\s-", "dung\\sNorway", "wood\\sHong", "dung, herbivore","unknown Unknown",
             "air\\s[A-Z][a-z]*", "soil,\\s[A-Z][a-z]*", "unknown,\\s[a-z]*", "unknown,\\s[A-z][a-z]*",
             "wood,\\s[a-z]*", "soil\\s[A-Z][a-z]*", "ground\\s[A-Z][a-z]*",
             "wood\\s[A-Z][a-z]*", "unknown\\s[A-Z][a-z]*", "water\\s[A-Z][a-z]*",
             "stem\\s[A-Z][a-z]*", "straw\\s[A-Z][a-z]*","paper", "leaves,\\s[A-Z][a-z]*",
             "paper\\sChina", "paper\\s[A-Z][a-z]*","wall\\spod", "rotten", "soil\\-",
             "\\s-","nodata", ":", "unknown\\s", "root\\s", "wood\\s", " f.\\ssp.", "BPI",
             "\\(White\\sspongy", "England", "Germany", "Finlland", "Ukraine", "United\\sStates",
             "Canada", "China", "Denmark", "hardwood", "conifer","deciduous", "conifer" )
  del <- paste(regex, collapse="|")
  x <- gsub(del, "", x)
  
  regex <- c("^\\s+", "\\s+$")
  del <- paste(regex, collapse="|")
  if(length(grep(del, x)) > 0) x <- gsub(del, "", x)

  x <- x[!(lapply(x, nchar) == 0 | lapply(x, nchar) == 1 | lapply(x, nchar) == 2)]
}
  return(x)
}
