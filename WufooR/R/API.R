#' @title Authentication
#' 
#' @description Method for setting your Wufoo Name and API Key. Your Wufoo name is the subdomain 
#' of your Wufoo URL displayed in "Account" tab. For example, for http://johnmalc.wufoo.com ->
#' the Wufoo_Name = "johnmalc". This can be also actually the company's name, e.g. \code{http://google.wufoo.com}
#' 
#' Your API may be found by selecting "Share" for your Form, then 
#' "API Information" (or go to \code{https://yourName.wufoo.com/api/code/1/}. 
#'  
#' @author The code for these methods has been developed by Scott Chamberlain \url{https://github.com/sckott} for his 
#' \url{https://github.com/ropensci/rnoaa} package. His copyright!
#' 
#' @param x - an empty parameter, e.g. NULL
#' 
#' @note Wufoo currently restricts your API usage to 5000 requests per day.
#'  
#' @examples 
#' options(Wufoo_Name = "johnmalc", Wufoo_API = "F1QH-Q64B-BSBI-JASJ")
#' 
#' @export
auth_name <- function(x) {
  tmp <- if(is.null(x)) {
    Sys.getenv("Wufoo_Name", "")
  } else x
  
  if(tmp == "") {
    getOption("Wufoo_Name", stop("you need to set up your wofoo name"))
  } else tmp
}

#' @rdname auth_name
#' @export
auth_key <- function(x) {
  tmp <- if(is.null(x)) {
    Sys.getenv("Wufoo_API", "")
  } else x
  
  if(tmp == "") {
    getOption("Wufoo_API",  stop("you need to set up your wofoo api key"))
  } else tmp
}

#' @title Generalized function for executing GET requests by always appending user's API Key.
#' 
#' @param url - which is used for the request
#' @param apiKey - uses the passed api key of the user
#' @param queryParameters - parameters that are used for building a URL
#' @param showURL - for debugging purposes only: it shows what URL has been called
#' @param debugConnection - same as above 
#' 
#' @import httr
#' @import jsonlite
#' 
#' @noRd
doRequest <- function(url, queryParameters = NULL, apiKey = auth_key(NULL), showURL = NULL, debugConnection = 0L) {
  
  if (is.null(apiKey)) {
    stop("Please assign your API Key", call. = FALSE)
  } else {
    
    # https://github.com/wufoo/Wufoo-PHP-API-Wrapper/blob/master/WufooApiWrapperBase.php#L102
    # http://curl.haxx.se/libcurl/c/CURLOPT_SSL_VERIFYHOST.html
    # https://stackoverflow.com/questions/28622558/how-to-solve-error-ssl23-get-server-hellosslv3-alert-handshake-failure

    if (.Platform$OS.type == "windows") {
      getResponse <- GET(url = url, query = queryParameters, 
                         config(userpwd = paste0(apiKey,":fakepassword"), ssl_cipher_list = "TLSv1", 
                                ssl_verifypeer=0L, ssl_verifyhost=0L, followlocation=1L, verbose=debugConnection))
    } else {
      getResponse <- GET(url = url, query = queryParameters, 
                         config(userpwd = paste0(apiKey,":fakepassword"), 
                                ssl_verifypeer=0L, ssl_verifyhost=0L, followlocation=1L, verbose=debugConnection))
      
    }
    stop_for_status(getResponse)
    
    rawTextResponse <- content(getResponse, as = "text")
    
    if (grepl("application/json", getResponse$headers$`content-type`)) {
      response <- fromJSON(rawTextResponse)
    } else {
      response <- rawTextResponse
    }
    
    if (identical(showURL, TRUE)) {
      cat("The requested URL has been this: ", getResponse$url, "\n") 
    }
    
    return(response)
  }
  
}
