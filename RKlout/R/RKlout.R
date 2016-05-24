## Source code for the package RKlout
## Author: Binayak Goswami
## Version: 1.0

## KLout API URLs
urlNetId <- "http://api.klout.com/v2/identity.json/twitter?screenName="
urlGetScore <- "http://api.klout.com/v2/user.json/"

## Function RKlout 
RKlout <- function(apiKey,twUser){
  apiGetNetId<-gsub("\\s","",paste(urlNetId,twUser,"&key=",apiKey))
  resultNetId <- RCurl::getURL(apiGetNetId)
  netId <- gsub("[^0-9.]","",resultNetId)
  apiGetKScore <- gsub("\\s","",paste(urlGetScore,netId,"/score?key=",apiKey))
  resultKScore <- strsplit(gsub("[^0-9.]"," ",RCurl::getURL(apiGetKScore)),split = " ")[[1]][10]
  nKloutScore <- as.numeric(resultKScore)
  ## Klout score can never be more than 100
  if(is.na(nKloutScore) || nKloutScore > 100 || is.null(nKloutScore))
  {
    errMsg <- "Please provide valid API key and Twitter Username/handle"
    return(errMsg)
  }
    return(nKloutScore)
}

