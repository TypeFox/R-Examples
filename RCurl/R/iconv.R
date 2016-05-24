
RCurlIconv =
function(str, from = "C99", to = "UTF-8", check = TRUE, quiet = FALSE)
{
  to = toupper(to)
  from = toupper(from)
  
  if(check) {
    w = c(from, to) %in% iconvlist()
    if(!all(w)) {
      if(!quiet)
        warning(paste(c(from, to)[w], collapse = ", "), " not supported iconv entries")
      return(str)
    }
  }
  ans = iconv(str, from, to)
  if(is.na(ans))
     str
  else
     ans
}
