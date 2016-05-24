#' Demo multipart parser with rhttpd
#'
#' Starts the Rhttpd web server and hosts a simple form including a file
#' upload to demo the multipart parser.
#
#' @export
#' @importFrom tools startDynamicHelp
demo_rhttpd <- function(){
  rhttpd_handler <- function(reqpath, reqquery, reqbody, reqheaders){

    # Extract HTTP content type and method from strange rhttpd format
    content_type <- grep("Content-Type:", strsplit(rawToChar(reqheaders), "\n")[[1]], ignore.case=TRUE, value=TRUE);
    content_type <- sub("Content-Type: ?", "", content_type, ignore.case=TRUE);
    http_method <- grep("Request-Method:", strsplit(rawToChar(reqheaders), "\n")[[1]], ignore.case=TRUE, value=TRUE);
    http_method <- sub("Request-Method: ?", "", http_method, ignore.case=TRUE);

    # Show HTML page for GET requests.
    if(http_method == "GET" || is.null(reqbody)){
      message("Received HTTP GET request: ", reqpath)
      testpage <- system.file("testpage.html", package="webutils");
      stopifnot(file.exists(testpage))
      list(
        "payload" = readBin(testpage, raw(), n=file.info(testpage)$size),
        "content-type" = "text/html",
        "headers" = NULL,
        "status code" = 200
      )
    } else {
      # Parse the multipart/form-data
      message("Received HTTP POST request.")

      # Check for multipart()
      postdata <- parse_http(reqbody, content_type)

      # Print it to the R console (just for fun)
      str(postdata)

      # process this form
      username <- rawToChar(as.raw(postdata$username$value))
      email <- rawToChar(as.raw(postdata$email_address$value))
      food <- rawToChar(as.raw(postdata$food$value))
      picture <- file.path(getwd(), basename(postdata$picture$filename))
      writeBin(as.raw(postdata$picture$value), picture)

      # return summary to the client
      list(
        "payload" = paste0("User: ", username, "\nEmail: ", email, "\nPicture (copy): ", picture,"\nFood: ", food, "\n"),
        "content-type" = "text/plain",
        "headers" = NULL,
        "status code" = 200
      )
    }
  }

  # Start rhttpd and get port
  port <- if(R.version[["svn rev"]] < 67550) {
    try(startDynamicHelp(TRUE), silent=TRUE);
    getFromNamespace("httpdPort", "tools");
  } else {
    startDynamicHelp(NA);
  }

  handlers_env <- getFromNamespace(".httpd.handlers.env", "tools")
  assign("test", rhttpd_handler, handlers_env)
  url <- paste0("http://localhost:", port, "/custom/test")
  message("Opening ", url)
  browseURL(url)
}
