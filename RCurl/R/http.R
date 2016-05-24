httpGET =
function(url, ..., curl = getCurlHandle())
{
  getURLContent(url, ..., curl = curl)
}

httpPOST =
function(url, ..., curl = getCurlHandle())
{
  getURLContent(url, .opts = list(customrequest = "POST", ...), curl = curl, post = 1L)
}

PUT = httpPUT =
function(url, content = NULL, ..., curl = getCurlHandle())
{
  if(!missing(content)) {
       val = if(is.character(content))
                charToRaw(paste(content, collapse = "\n"))
             else if(is.raw(content))
                content
             else
                stop("not certain how to convert content to the target type for a PUT request")
       
       getURLContent(url, infilesize = length(val),
                          readfunction = val,
                          upload = TRUE, 
                          ..., curl = curl, customrequest = "PUT")
  } else
      getURLContent(url, ..., curl = curl, customrequest = "PUT")
}

DELETE = httpDELETE =
function(url, ..., curl = getCurlHandle())
{
  getURLContent(url, customrequest = "DELETE", ..., curl = curl)
}

HEAD = httpHEAD =
function(url, ..., curl = getCurlHandle())
{
  getURLContent(url, customrequest = "HEAD", nobody = TRUE, ..., curl = curl)
}

httpOPTIONS =
function(url, ..., curl = getCurlHandle())
{
  ans = getURLContent(url, customrequest = "OPTIONS", ..., curl = curl, header = TRUE)
  ans$header
}


