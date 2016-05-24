getUserName =
function()
{
   Sys.getenv(if(.Platform$OS.type == "windows") "USERNAME" else "USER")
}

scp =
function(host, path, keypasswd = NA, user = getUserName(), rsa = TRUE,
         key = sprintf(c("~/.ssh/id_%s.pub", "~/.ssh/id_%s"), if(rsa) "rsa" else "dsa"),
         binary = NA, size = 5000, curl = getCurlHandle(), ...)
{
  if(!is.na(user) && length(grep("@", host)) == 0)
    host = paste(user, "@", host, sep = "")
  
  if(length(grep("^scp", host)) == 0)
    host = paste("scp://", host, sep = "")

  if(!missing(path))
    host = paste(host, path, sep = "/")

  key = path.expand(key)
  opts = list(ssh.public.keyfile = key[1], ssh.private.keyfile = key[2], ...)
  if(!is.na(keypasswd))
    opts$keypasswd = keypasswd

  # We don't use getURLContent() as it wants to process the HTTP header
  # which isn't part of this. We could provide our own version of the header
  # argument but this is essentially what we are doing by specifying writefunction.
           #  getURLContent(host, .opts = opts)
  opts$url = host

  if(is.na(binary) || binary) {
    binary = TRUE
    opts$writefunction = getNativeSymbolInfo("R_curl_write_binary_data")$address
    buf <- binaryBuffer(size)    
    opts$file = buf@ref
  } else {
    h = basicTextGatherer()
    opts$writefunction = h$update
  }
  curlPerform(.opts = opts, curl = curl)

  if(binary)
     as(buf, "raw")
  else
     h$value()
}
