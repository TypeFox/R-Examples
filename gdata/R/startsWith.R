startsWith <- function(str, pattern, trim=FALSE, ignore.case=FALSE)
  {
    if(trim) str <- trim(str)
    if(ignore.case)
      {
        str <- toupper(str)
        pattern <- toupper(pattern)
      }
    substr(str,start=1,stop=nchar(pattern))==pattern
  }
