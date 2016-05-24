setRefClass("OAuth",
            fields = list(
              consumerKey = "character",
              consumerSecret = "character",
              oauthKey = "character",
              oauthSecret = "character",
              needsVerifier = "logical",
              handshakeComplete = "logical",
              verifier = "character",
              requestURL = "character",
              authURL = "character",
              accessURL = "character",
              signMethod = 'character'
              ),
            methods = list(
              initialize = function(needsVerifier, ...) {
                if (!missing(needsVerifier))
                  needsVerifier <<- needsVerifier
                else
                  needsVerifier <<- TRUE
                handshakeComplete <<- FALSE
                callSuper(...)
                .self
              },
              
              handshake = function(signMethod='HMAC', curl=getCurlHandle(), browseUrl=TRUE, ...) {
                ' Performs the OAuth handshake.  In most cases
                  the user will need to complete a manual step
                  with their web browser, entering a PIN into
                  this function.
                '
                .self$handshakeComplete <- FALSE
                .self$signMethod <- signMethod
                resp <- oauthPOST(.self$requestURL, 
                                  .self$consumerKey,
                                  .self$consumerSecret,
                                  NULL, NULL, signMethod=.self$signMethod, curl=curl,
                                  handshakeComplete=.self$handshakeComplete, ...)
                vals <- parseResponse(resp)
                if (!all(c('oauth_token', 'oauth_token_secret') %in%
                         names(vals))) {
                  stop("Invalid response from site, please ",
                     "check your consumerKey and consumerSecret",
                     " and try again.")
                }
                .self$oauthKey <- vals["oauth_token"]
                .self$oauthSecret <- vals["oauth_token_secret"]
                if (.self$needsVerifier) {
                  verifyURL <- paste(.self$authURL, "?oauth_token=",
                                     oauthKey, sep='')
                  msg <- paste("To enable the connection, please direct",
                               " your web browser to: \n",
                               verifyURL,
                               "\nWhen complete, record the PIN given ",
                               "to you and provide it here: ", sep='')
                  if (browseUrl) {
                    browseURL(verifyURL)
                  }
                  .self$verifier <- readline(prompt=msg)
                }
                params <- c(oauth_verifier=.self$verifier)
                resp <- oauthPOST(.self$accessURL, .self$consumerKey, .self$consumerSecret,
                                  .self$oauthKey, .self$oauthSecret, signMethod=.self$signMethod,
                                  curl=getCurlHandle(), params=params, handshakeComplete=.self$handshakeComplete,
                                  ...)
                vals <- parseResponse(resp)
                if (!all(c('oauth_token', 'oauth_token_secret') %in%
                         names(vals))) {
                  stop("Invalid response after authorization.  ",
                       "You likely misentered your PIN, try rerunning",
                       " this handshake & browser authorization to get",
                       " a new PIN.")
                }
                .self$oauthKey <- vals["oauth_token"]
                .self$oauthSecret <- vals["oauth_token_secret"]
                .self$handshakeComplete <- TRUE
              },
              
              isVerified = function() {
                'Will report if this object is verified or not.
                 Verification can either involve not needing it
                 in the first place, or as part of the handshake'
                if (.self$needsVerifier)
                  length(verifier) != 0
                else
                  TRUE
              },
              
              OAuthRequest = function(URL, params=character(), method="GET",
                                      customHeader=NULL, curl=getCurlHandle(), ...) {
                ' If the OAuth handshake has been completed, will
                submit a URL request with an OAuth signature, returning
                any response from the server
                '
                if (! .self$handshakeComplete)
                  stop("This OAuth instance has not been verified")

                httpFunc <- switch(method,
                                   POST = oauthPOST,
                                   GET = oauthGET,
                                   stop("method must be POST or GET"))

                ## RCurl converts `\\` to `\` if the content includes strings like `\u0`
                ## and for example, response JSON is broken as a result,
                ## so change RCurl:::mapUnicodeEscapes temporarily to avoid it
                ## cf. https://github.com/omegahat/RCurl/issues/1
                ##rcurlEnv <- getNamespace("RCurl")
                ##mapUnicodeEscapes <- get("mapUnicodeEscapes", rcurlEnv)
                ##unlockBinding("mapUnicodeEscapes", rcurlEnv)
                ##assign("mapUnicodeEscapes", function(str) str, rcurlEnv)
                ##on.exit({
                ##  assign("mapUnicodeEscapes", mapUnicodeEscapes, rcurlEnv)
                ##  lockBinding("mapUnicodeEscapes", rcurlEnv)
                ##}, add = TRUE)

                httpFunc(URLencode(URL), params=params, consumerKey=.self$consumerKey,
                         consumerSecret=.self$consumerSecret,
                         oauthKey=.self$oauthKey, oauthSecret=.self$oauthSecret,
                         customHeader.self$customHeader, curl=getCurlHandle(), signMethod=.self$signMethod, ...)
              }
              )
            )

OAuthFactory <- getRefClass("OAuth")
OAuthFactory$accessors(names(OAuthFactory$fields()))


parseResponse <- function(response) {
  ## Will return a named vector, so a response field of the
  ## form foo=blah&qwerty=asdf will have vals c(blah,asdf) and
  ## names will be c(foo, qwerty).  If the response is borked,
  ## the output of this function will be as well, so caveat
  ## emptor, GIGO, etc
  pairs <- sapply(strsplit(response, '&')[[1]], strsplit, '=')
  out <- sapply(pairs, function(x) x[2])
  names(out) <- sapply(pairs, function(x) x[1])
  out
}

oauthPOST <- function(url, consumerKey, consumerSecret,
                      oauthKey, oauthSecret, params=character(), customHeader = NULL,
                      curl = getCurlHandle(), signMethod='HMAC', handshakeComplete=TRUE, ...) {
  if(is.null(curl))
    curl <- getCurlHandle()
  
  params <- signRequest(url, params, consumerKey, consumerSecret,
                      oauthKey=oauthKey, oauthSecret=oauthSecret,
                      httpMethod="POST", signMethod=signMethod,
                      handshakeComplete=handshakeComplete)
  opts <- list(...)
  
  ## post ,specify the method
  postForm(url, .params = params, curl = curl,
           .opts = opts, style = "POST")
}

#XXX? use .opts for the curl options.
#     add a ... for the parameters.

oauthGET <- function(url, consumerKey, consumerSecret,
                     oauthKey, oauthSecret, params=character(), customHeader = NULL,
                     curl = getCurlHandle(), signMethod='HMAC', ..., .opts = list(...)) {
  ##   opts = list(httpheader = c(Authentication = paste(names(auth),  auth, sep = "=", collapse = "\n   ")), ...)
  if(is.null(curl))
    curl <- getCurlHandle()
   
   params <- signRequest(url, params, consumerKey, consumerSecret,
                       oauthKey=oauthKey, oauthSecret=oauthSecret,
                       httpMethod="GET", signMethod=signMethod)

   getForm(url, .params = params, curl = curl, .opts = c(httpget = TRUE,  list(...)))
}
