# This function is for parsing the parameters in a URL
# for a form.
# It is based on mail and a function from Chris Davis (Delft Uni.)
#

getFormParams =
  #
  #  getFormParams("http://www.omegahat.net/foo/bob.R?xyz=1&abc=verylong")
  #  getFormParams("xyz=1&abc=verylong")
  #  getFormParams("xyz=1&abc=")
  #  getFormParams("xyz=1&abc=&on=true")
  #  getFormParams("")
  #  getFormParams(character())
   # This handles but doesn't detect that there is no =
  #  getFormParams("abc")
  
function(query, isURL = grepl("^(http|\\?)", query))
{
  if(length(query) == 0)
    return(NULL)
  
       # allow for either a full URL with the query parameters
       # or alternative, just the query parameters.
  if(isURL)     # use parseURI from the XML package, but could
                # use gsub().
                # parseURI(query)$query
     query = gsub(".*\\?", "", query)

  if(nchar(query) == 0)
    return(character())

  els = strsplit(query, "[&=]")[[1]]
  i = seq(1, by = 2, length = length(els)/2)
  ans = structure(els[i+1L], names = els[i])
    # if the last parameter didn't have a value, we'll get an NA
    # Not any of the the internal ones.
  if(any(i <- is.na(ans)))
     ans[i] = ""

  names(ans) = trim(names(ans))
  ans
}

trim =
function(x)
  gsub("(^[[:space:]]+|[[:space:]]$)", "", x)


