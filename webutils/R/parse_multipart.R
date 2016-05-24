#' Parse a multipart/form-data request
#'
#' Parse a multipart/form-data request, which is usually generated from a HTML form
#' submission. The parameters can include both text values as well as binary files.
#' They can be distinguished from the presence of a \code{filename} attribute.
#'
#' A multipart/form-data request consists of a single body which contains one or more
#' values plus meta-data, separated using a boundary string. This boundary string
#' is chosen by the client (e.g. the browser) and specified in the \code{Content-Type}
#' header of the HTTP request. There is no escaping; it is up to the client to choose
#' a boundary string that does not appear in one of the values.
#'
#' The parser is written in pure R, but still pretty fast because it uses the regex
#' engine.
#'
#' @export
#' @param body body of the HTTP request. Must be raw or character vector.
#' @param boundary boundary string as specified in the \code{Content-Type} request header.
#' @examples \dontrun{example form
#' demo_rhttpd()
#' }
parse_multipart <- function(body, boundary){
  # Some HTTP daemons give the body as a string instead of raw.
  if(is.character(body))
    body <- charToRaw(paste(body, collapse=""))

  if(is.raw(boundary))
    boundary <- rawToChar(boundary)

  stopifnot(is.raw(body))
  stopifnot(is.character(boundary))

  # Extract the boundary string from the HTTP header
  boundary <- paste0("--", boundary)
  boundary_length <- nchar(boundary);

  # Find the locations of the boundary string
  indexes <- grepRaw(boundary, body, fixed=TRUE, all=TRUE);

  if(!length(indexes))
    stop("Boundary was not found in the body.")

  if(length(indexes) == 1){
    if(length(body) < (boundary_length+5)){
      # Empty HTML5 FormData object
      return(list())
    } else {
      # Something went wrong
      stop("The 'boundary' was only found once in the multipart/form-data message. It should appear at least twice. The request-body might be truncated.")
    }
  }

  parts <- list();
  for(i in seq_along(head(indexes, -1))){
    from <- indexes[i] + boundary_length;
    to <- indexes[i+1] -1;
    parts[[i]] <- body[from:to];
  }

  out <- lapply(parts, multipart_sub);
  names(out) <- vapply(out, function(x){as.character(x$name)}, character(1))
  return(out)
}

multipart_sub <- function(bodydata){
  stopifnot(is.raw(bodydata));
  splitchar <- grepRaw("\\r\\n\\r\\n|\\n\\n|\\r\\r", bodydata);
  if(!length(splitchar))
    stop("Invalid multipart subpart:\n\n", rawToChar(bodydata));

  headers <- bodydata[1:(splitchar-1)];
  headers <- trail(rawToChar(headers));
  headers <- gsub("\r\n", "\n", headers);
  headers <- gsub("\r", "\n", headers);
  headerlist <- unlist(lapply(strsplit(headers, "\n")[[1]], trail));

  dispindex <- grep("^Content-Disposition:", headerlist);
  if(!length(dispindex))
    stop("Content-Disposition header not found:", headers);
  dispheader <- headerlist[dispindex];

  #get parameter name
  m <- regexpr("; name=\\\"(.*?)\\\"", dispheader);
  if(m < 0)
    stop('failed to find the name="..." header')

  namefield <- unquote(sub("; name=", "", regmatches(dispheader, m), fixed=TRUE));

  #test for file upload
  m <- regexpr("; filename=\\\"(.*?)\\\"", dispheader);
  if(m < 0){
    filenamefield = NULL;
  } else {
    filenamefield <- unquote(sub("; filename=", "", regmatches(dispheader, m), fixed=TRUE));
  }

  #filedata
  splitval  <- grepRaw("\\r\\n\\r\\n|\\n\\n|\\r\\r", bodydata, value=TRUE);
  start <- splitchar + length(splitval);
  if(identical(tail(bodydata,2), charToRaw("\r\n"))){
    end <- length(bodydata)-2;
  } else {
    end <- length(bodydata)-1;
  }

  #the actual fields
  list (
    name = namefield,
    value = bodydata[start:end],
    filename = filenamefield
  )
}
