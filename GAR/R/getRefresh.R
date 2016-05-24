getRefresh <- 

function(client_id=Sys.getenv('GAR_CLIENT_ID'), client_secret=Sys.getenv('GAR_CLIENT_SECRET'), code=code){
  
        ###POST TO GET REFRESH TOKEN
        rToken <- POST(
                      "https://accounts.google.com/o/oauth2/token",
                      body = list(
                                  "code"=code,
                                  "client_id"=client_id,
                                  "client_secret"=client_secret,
                                  "redirect_uri"="http://localhost",
                                  "grant_type"="authorization_code"
                                  ),
                      add_headers('Host'='accounts.google.com'), 
                      encode='form'
                      )
        
        ##PARSE REFRESH TOKEN RESPONSE
        rToken <- content(rToken)
        rToken <- rToken$refresh_token
        return(rToken)

}