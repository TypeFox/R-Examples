# starts session on avis server

.ravis_session_started<- list( OK = 1, NO_COOKIES = 2, BAD_CREDENTIALS = 3 )

.ravis_url_login <- "http://proyectoavis.com/cgi-bin/login.cgi"

# Get contents of a url. Depending on 'nologin' parameter it uses an avisproject.com user session
# or not
.avisGetURL <- function(url, nologin = FALSE) {
  if (nologin == TRUE){
    # new curl handle
    curl_handler<- getCurlHandle()
  } else {
    curl_handler<- .avisCurlHandler()
  }

  return (getURL(url, curl = curl_handler))
}

.avisCurlHandler <- function(){

  if(!.avisCacheHas(".ravis_curl_handler")){
    .avisUserLogin()
  }

  return (.avisCacheGet(".ravis_curl_handler"))
}

# logs to web with ravis specific user
.avisUserLogin <- function() {

    # hack for avoiding NOTE on check: 'no visible binding for global variable'
    # see: http://stackoverflow.com/questions/9439256/how-can-i-handle-r-cmd-check-no-visible-binding-for-global-variable-notes-when
    ravis_credentials <- NULL
    rm(ravis_credentials)
    
    return (.avisLogin("ravis-user", ravis_credentials[[1]]))
}

# logs to web with user defined credentials
.avisLogin <- function (avis_user, avis_pass) {
  # log user in remote server

  params<- list( usu=avis_user, password=avis_pass, control_login='1' )

  html = postForm(.ravis_url_login, 
    .params = params, 
    curl = .avisCreateCurlHandler(), 
    style="POST",
    .encoding="ISO-8859-1")

  status<- .parse.avisLoginStatusFromHTML(html)[1]

  if(status == .ravis_session_started["NO_COOKIES"]){
    stop("There is no cookie support in the request. Unable to initializate session")
  } else if(status == .ravis_session_started["BAD_CREDENTIALS"]) {
    stop(paste("Wrong username and/or password for user: ", avis_user))
  }

  return (status == .ravis_session_started["OK"])
}

# Create new Curl handle for the requests to avis
.avisCreateCurlHandler <- function() {
  
  .avisVerboseMessage("INFO: initializing curl handler for connecting to avis project")

  .handler <- getCurlHandle()

  curlSetOpt(
    .opts = list(cainfo = system.file("CurlSSL", "cacert.pem", package = "RCurl")),
      cookiefile = tempfile("r_avis_cookie.txt"),
      useragent = 'R-Avis package',
      followlocation = TRUE,
      httpheader = "Referer: http://proyectoavis.com",
      curl = .handler)

  .avisCacheSet(".ravis_curl_handler", .handler)

  # first call to initializate session
  getURL(.ravis_url_login, curl = .handler)

  return (.handler)
}

# Search HTML to find out the result of the login request
.parse.avisLoginStatusFromHTML<- function(html) {
  
  html <- tolower(html)

  status <- NULL
  
  if(.textHasString(html, "usuario o clave incorrecta")){
    status <- .ravis_session_started["BAD_CREDENTIALS"]
  }

  if(.textHasString(html, "navegador no acepta cookies")){
    status <- .ravis_session_started["NO_COOKIES"]
  }

  if(.textHasString(html, "bienvenido")){
    status <- .ravis_session_started["OK"]
  }

  if(is.null(status)){
    stop("Could not find out the current login status")
  }

  return (as.integer(status[[1]]))
}

# Checks whether a text contains some string, returning logical
.textHasString <- function(text, string) {
  
  res <- grep(string, text)

  if(1 == length(res) && 1 == res[1]){
    return (TRUE)
  }

  return (FALSE)
}