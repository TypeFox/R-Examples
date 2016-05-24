mixpanelCreateAccount <- function(
  name,         # Arbitrary name.
  token,        # Mixpanel token.
  key,          # API key of the Mixpanel project.
  secret,       # API Secret of the Mixpanel project.
  mongoDBname,  # Used when imported into a MongoDB database. 
  dataPath,     # File path to store exported raw data.
  RDataPath     # File path to store R raw data.
  ) {
  obj = list(name=name, token=token, apiKey=key, apiSecret=secret)
  
  if (!missing(mongoDBname))
    obj$mongoDBname=mongoDBname
  
  if (!missing(dataPath))
    obj$dataPath = dataPath
  
  if (!missing(RDataPath))
    obj$RDataPath = RDataPath
  
  class(obj) = "mixpanelAccount"
  obj
}
