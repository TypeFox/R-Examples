#' A class for managing Google OAuth2 workflow for installed apps
#' 
#' This class is aimed to manage the process where the user
#' grants Google account access to your R application.
#'
#' There are two basic ways to create an GoogleOAuth2 object.
#' First, you can use \code{google.oauth2} to initially register
#' your application with google. You may store this authentication
#' information using \code{put}. Second, you can use \code{read.google.oauth2}
#' to load authentication information from a previously written file.
#'
#' @seealso \code{\link{google.oauth2}}, 
#'
#' @references 
#' \href{https://developers.google.com/accounts/docs/OAuth2InstalledApp?hl=de}{Google}
#' @name GoogleOAuth2-class
#' @rdname GoogleOAuth2-class
#' @exportClass GoogleOAuth2
setClass(
    Class="GoogleOAuth2", 
    representation=representation(
        access_token="character",
        access_expires="POSIXct",
        refresh_token="character",
        scope="character",
        client_id="character", 
        client_secret="character"
    ), 
    contains="Target"
)


#' Google OAuth2 initial authentication
#'
#' This function redirects the user to a consent screen where the user can
#' grant or revoke access for the provided scope.
#'
#' @param scope          either shortname ("blogger") or URL of the scope
#' @param client_id      client id. Defaults to getOption("datamart.client_id")
#' @param client_secret  client secret. Defaults to getOption("datamart.client_secret")
#' @param name           currently ignored
#' @param clss           the class to create. Must be derived from GoogleOAuth2
#' @param curl           curl.handle object
#' @param verbose        prints diagnostic messages on the way
#'
#' @seealso \code{\link{read.google.oauth2}}, \code{\link{google.oauth2}}
#'
#' @return GoogleOAuth2 object 
#' @export
google.oauth2 <- function(
    scope, 
    client_id=getOption("datamart.client_id"), 
    client_secret=getOption("datamart.client_secret"), 
    name="", 
    clss="GoogleOAuth2",
    curl=RCurl::getCurlHandle(),
    verbose=TRUE
) {
    uri <- 'https://accounts.google.com/o/oauth2/auth'
    scopes <- c(blogger="https://www.googleapis.com/auth/blogger")
    
    if(!grepl('^https', scope))  {
        if(scope %in% names(scopes)) 
            scope <- scopes[[tolower(scope)]]
        else
            stop("invalid scope for google oauth2.")
    }

    if(is.null(client_id)) stop("missing client_id for google oauth2.")
    if(is.null(client_secret)) stop("missing client_secret for google oauth2.")
    
    args <- c(
        redirect_uri = "urn:ietf:wg:oauth:2.0:oob",
        scope = scope,
        client_id = client_id,
        response_type = 'code'
    )
    uri = sprintf("%s?%s", uri, paste(names(args), args, sep = "=", collapse = '&'))
    browseURL(uri)
    cat("Cut-and-paste the permission string from your Web browser here: ")
    code <- readLines(stdin(), 1)
    
    ret <- RCurl::postForm(
        'https://accounts.google.com/o/oauth2/token',
        .params=list(
            code=code,
            client_id=client_id,
            client_secret=client_secret,
            redirect_uri = "urn:ietf:wg:oauth:2.0:oob",
            grant_type="authorization_code"
        ),
        style="POST",
        header=TRUE,
        curl=curl,
        binary=FALSE,
        .opts=list(
            # followlocation=TRUE,
            ssl.verifypeer = TRUE, 
            ssl.verifyhost = TRUE, 
            cainfo = system.file("CurlSSL", "cacert.pem", package = "RCurl")
        ),
        .contentEncodeFun=RCurl::curlPercentEncode
    )
    
    if(!attr(ret, "Content-Type")=="application/json") stop("initial google oauth2 authentication failed.")
    
    tokens <- RJSONIO::fromJSON(ret)
    if(verbose) cat("Access-Token: ", tokens[["access_token"]],"\n")

    new(
        clss, 
        name=name, 
        access_token=tokens[["access_token"]],
        access_expires=Sys.time()+tokens[["expires_in"]]-1,
        refresh_token=tokens[["refresh_token"]],
        scope=scope,
        client_id=client_id, 
        client_secret=client_secret
    )
}

#' Refreshed Google OAuth2 authentication
#'
#' This function reads a file previously created via \code{put(oauth2, filename)}.
#'
#' @param fn   data filename
#'
#' @seealso \code{\link{google.oauth2}}
#'
#' @return GoogleOAuth2 object 
#' @export
read.google.oauth2 <- function(fn) {
    e <- new.env()
    load(fn, envir=e)
    e[["target"]]
}

#' @param curl    (GoogleOAuth2) curl handle
#' @rdname query-methods
#' @name query
#' @export
#' @docType methods
#' @aliases query query,GoogleOAuth2,character-method
setMethod(
    f="query",
    signature=c(self="GoogleOAuth2", resource="character"),
    definition=function(self, resource, curl=RCurl::getCurlHandle(), ...) {
        if(resource=="access_token") {
            if(self@access_expires<Sys.time()) {
                ret <- RCurl::postForm(
                    'https://accounts.google.com/o/oauth2/token',
                    .params=list(
                        refresh_token=self@refresh_token,
                        client_id=self@client_id,
                        client_secret=self@client_secret,
                        grant_type="refresh_token"
                    ),
                    style="POST",
                    header=TRUE,
                    curl=curl,
                    binary=FALSE,
                    .opts=list(
                        ssl.verifypeer = TRUE, 
                        ssl.verifyhost = TRUE, 
                        cainfo = system.file("CurlSSL", "cacert.pem", package = "RCurl")
                    ),
                    .contentEncodeFun=RCurl::curlPercentEncode
                )
    
                if(!attr(ret, "Content-Type")=="application/json") stop("refresh google oauth2 authentication failed.")
                tokens <- RJSONIO::fromJSON(ret)
                self@access_token <- tokens[["access_token"]]
                self@access_expires <- Sys.time()+tokens[["expires_in"]]-1
            }
            return(self@access_token)
        } else callNextMethod()
    }
)

#' @rdname put-methods
#' @name put
#' @export
#' @docType methods
#' @aliases put put,GoogleOAuth2,character-method
setMethod(
  f="put",
  signature=c(target="GoogleOAuth2", where="character"),
  definition=function(target, where, ...) save(target, file=where)
)
