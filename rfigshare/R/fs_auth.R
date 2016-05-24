
#' Figshare authentication via OAuth1.0 using httr
#'
#' Figshare authentication via OAuth1.0 using httr
#' @param cKey optional argument for the consumer key.  See details.  
#' @param cSecret optional argument for the consumer secret key. See details. 
#' @param token optional argument for the user-specific token.  See details.  
#' @param token_secret Optional argument to provide a secret token assigned
#' to the user, rather than let fs_auth() automatically handle authentication. 
#' See details.  
#' @details Explicit calls to fs_auth() are usually not needed,
#' as the function is called automatically by all other functions that
#' need authentication.  As of version 0.3, no arguments are needed as
#' authentication is done online, and fs_auth() will not attempt to load
#' keys stored in options.  
#' 
#' By default, the function will use the application's consumer key and
#' consumer secret key, rather than expecting the user to create their own
#' application.  The user-specific tokens will then be generated and locally
#' cached for use between sessions, if indicated by the interactive options. 
#' For details, see httr oauth1.0_token documentation.  
#'
#' If for some reason a user would rather provide there token and secret token
#' as before this is still supported using the same arguments.  Users wanting to
#' have their own app can provide cKey and cSecret arguments too, but this
#' is provided primarily for backwards compatibility with older versions.  It
#' is expected that most users will leave the keys as NULL.  
#' @return OAuth credential (invisibly).  The credential is also written to the enivronment "FigshareAuthCache", which is created when the package is loaded.  All functions needing authentication can then access the credential without it being explicitly passed in the function call. If authentication fails, returns the failing GET response for debugging.   
#' @import httr httpuv
#' @export
#' @author Carl Boettiger \email{cboettig@@gmail.com}
#' @references \url{http://api.figshare.com}
#' @examples \dontrun{
#' fs_auth()
#' }
fs_auth <- 
  function(cKey = getOption("FigshareKey", NULL),
           cSecret = getOption("FigsharePrivateKey", NULL),
           token = getOption("FigshareToken", NULL),
           token_secret = getOption("FigsharePrivateToken", NULL)){

  if(is.null(cKey))
    cKey <- "Kazwg91wCdBB9ggypFVVJg"
  if(is.null(cSecret))
    cSecret <- "izgO06p1ymfgZTsdsZQbcA"

  endpoint <- 
    oauth_endpoint("request_token", 
                   "authorize", 
                   "access_token", 
                   base_url = "http://api.figshare.com/v1/pbl/oauth")
  myapp <- oauth_app("rfigshare", 
                     key = cKey, 
                     secret = cSecret)

  if(is.null(token) && is.null(token_secret)) {
    oauth <- oauth1.0_token(endpoint, myapp)
  } else {
    resp <- sign_oauth1.0(myapp, token = token, token_secret = token_secret)
    oauth <- resp$token
  }

  assign('oauth', oauth, envir=FigshareAuthCache)

  # Test that we have authenticated
  response <- GET("http://api.figshare.com/v1/my_data/articles", 
         config(token = oauth)
         )

  if(response$status_code != 200){
    warning(paste("Authentication failed, please check your credentials. Got code",
                  response$status_code, "with response", content(response)))
    oauth <- response
  } else {
    message("Authentication successful")
  }
  invisible(oauth)
} 


