tokenRefresh <-
function(client_id=Sys.getenv('GAR_CLIENT_ID'), client_secret=Sys.getenv('GAR_CLIENT_SECRET'), token=Sys.getenv('GAR_REFRESH_TOKEN')) {
                response <- POST(
                                "https://accounts.google.com/o/oauth2/token",
                                body = list(
                                       "client_id"=client_id,
                                       "client_secret"=client_secret,
                                       "refresh_token"=token,
                                       "grant_type"="refresh_token"),
                                add_headers('Host'='accounts.google.com'), 
                                encode='form'
                                )
                response <- content(response)
                accessToken <- response$access_token
                assign('GAR_ACCESS_TOKEN', accessToken, envir=envGAR)
                return(accessToken)
                }


