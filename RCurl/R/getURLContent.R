getURLContent =
  #
  # Used to be
  #  header = basicTextGatherer()
  #  ans = getBinaryURL(url, headerfunction = header$update, curl = header$curl())
  #  processContent(ans, header$header(), .encoding)
  # but now we use the dynamic reader.
  #
function(url, ..., curl = getCurlHandle(.opts = .opts), .encoding = NA, binary = NA, .opts = list(...),
         header = dynCurlReader(curl, binary = binary, baseURL = url, isHTTP = isHTTP, encoding = .encoding),
          isHTTP = length(grep('^[[:space:]]*http', url)) > 0)
{
  url = as(url, "character")

  if(!missing(curl))
     curlSetOpt(.opts = .opts, curl = curl)

  if(is.logical(header)) {
     returnHeader = header
     header = dynCurlReader(curl, binary = binary, baseURL = url, isHTTP = isHTTP, encoding = .encoding)
  } else
     returnHeader = FALSE

  if(!('headerfunction' %in% names(.opts))) {
    # .opts$headerfunction = header$update
     protect = missing(header)
     curlSetOpt(curl = curl, .isProtected = protect,
                 headerfunction = header$update)
   }

  if(!isHTTP && !('writefunction' %in% names(.opts))) {
      # If for example this is scp where there is no header
      # or headerfunction will never get called. So we have to
      # set the writefunction as well.
    # .opts$headerfunction = header$update
     protect = missing(header)
     curlSetOpt(curl = curl, .isProtected = protect,
                 writefunction = header$update)
   }

  curlPerform(url = url, curl = curl, .opts = .opts)

  if(isHTTP && length(header$header())) {
     http.header = parseHTTPHeader(header$header())
     stop.if.HTTP.error(http.header)
  }

  if(returnHeader)
    list(header = if(is(returnHeader, "AsIs"))
                      header$header()
                  else
                      parseHTTPHeader(header$header()),
         body = header$value())
  else
    header$value()
}

stop.if.HTTP.error =
function(http.header)
{

  if(length(http.header) == 0)
    return(NA) # or TRUE

  if( floor(as.integer(http.header[["status"]])/100) >= 4) {
     klass = getHTTPErrorClass(http.header[["status"]])
     err = simpleError(http.header[["statusMessage"]])
     err$httpHeader = http.header
     class(err) = c(klass, class(err))
     #signalCondition(err)
     stop(err)
  }

  TRUE
}

processContent =
#
# Figure out how to interpret the contents based on the HTTP response's header
# i.e. look at its Content-Type.
#
function(ans, header, .encoding = NA)
{
  headerText = if(is.character(header)) header else header$value()
  http.header = parseHTTPHeader(headerText)

  stop.if.HTTP.error(http.header)

  content.type = getContentType(http.header)
  binary = isBinaryContent(http.header, content.type)
  if(!(is.na(binary) || binary)) {
     ans = rawToChar(ans)
     if(length(.encoding)  == 0 || is.na(.encoding)) {
        charset = grep("charset", content.type, value = TRUE)
        if(length(charset))
          .encoding = strsplit(charset, "=")[[1]][2]
     }
     if(length(.encoding) && !is.na(.encoding))
       Encoding(ans) = .encoding
  } else {
     attr(ans, "Content-Type") = getContentType(http.header)
     ans
  }
  ans
}

trim =
function(x)
{
    gsub("(^[[:space:]]+|[[:space:]]+$)", "", x, perl = TRUE)
}

getContentType =
function(header, full = FALSE)
{
   i = match("content-type", tolower(names(header)))
   if( is.na( i ) )
      return(character())

   tmp = trim(strsplit(header[i[1]], "; *")[[1]])
   if(!full)
       return(tmp)

   vals = strsplit(tmp, "=")
   structure(gsub(";$", "", sapply(vals, function(x) x[length(x)])),
             names = sapply(vals, function(x) if(length(x) > 1) x[1] else ""))
}

# See http://www.iana.org/assignments/media-types/
textContentTypes = c("html", "text", "xhtml", "plain", "xml", "x-latex", "css", "latex",
                     "sgml", "postscript", "texinfo", "ecmascript", "javascript",
                     "atom+xml", "json")

isBinaryContent =
  #
  #  type can be given as a list intended to be separate  header elements
  #  e.g. Content-Type, Content-Encoding, etc.
  #  Each can be a vector.
  #
function(header, type = getContentType(header)[1],
          textTypes = getOption("text.content.types"))
{
   if(is.list(type) && length(type) > 1) {
     last <- TRUE
     for(i in type) {
        if(length(i) && !is.na(i) && (last <- isBinaryContent(header, i, textTypes)))
          return(TRUE)
     }
     return(last)
   }

   if(length(type) == 0)
     return(NA)

   if(is.null(textTypes))
     textTypes = textContentTypes
   type.els = strsplit(type, "/")[[1]]
   if(type.els[1] == "text")
     return(FALSE)

   if(any(type.els %in% textContentTypes))
      return(FALSE)

   if(length(grep("\\+xml$", type.els)))
      return(FALSE)

   TRUE
}

