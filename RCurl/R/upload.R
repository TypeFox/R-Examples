ftpUpload =
  #
  # what is the name of the file or the contents of the file
  #
function(what, to, asText = inherits(what, "AsIs") || is.raw(what),
          ..., curl = getCurlHandle())
{
  if(!asText && !inherits(what, "connection")) {
    file = file(what, "rb")
    on.exit(close(file))
  } else
     file = what
  
  curlPerform(url = to, upload = TRUE,
              readfunction = uploadFunctionHandler(file, asText), ..., curl = curl)
}

uploadFunctionHandler =
  #
  # returns the function that is called as the READFUNCTION callback.
  #
  # This handles raw, character and file contents.
  #
function(file, asText = inherits(file, "AsIs") || is.raw(file))
{

  if(asText) {
      pos = 1
      isRaw = is.raw(file)
      len = if(isRaw) length(file) else nchar(file)

      function(size) {

        if(pos > len)
          return(if(isRaw) raw(0) else character(0))

        ans = if(isRaw)
                  file[seq(pos, length = size)]
               else
                  substring(file, pos, pos + size)
        pos <<- pos + len
        ans
      }

   } else

      function(size) {
         readBin(file, raw(), size)
      }
}


CFILE = 
function(filename, mode = "r")
{
   filename = path.expand(filename)
   .Call("R_openFile", filename, as.character(mode))
}


setMethod("close", "CFILE",
           function(con, ...) {
             .Call("R_closeCFILE", con@ref, PACKAGE = "RCurl")
             con
          })

