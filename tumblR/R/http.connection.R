http.connection <-
function(url,token,bodyParams,consumer_key,consumer_secret,connection="GET"){
  
  urlParsed <- parse_url(url)
  base_url <- build_url(urlParsed)
  
  oauth <- compact(list(
    oauth_consumer_key = token$app$key,
    oauth_nonce = nonce(),
    oauth_signature_method = "HMAC-SHA1",
    oauth_timestamp = as.integer(Sys.time()),
    oauth_version = "1.0",
    oauth_token = token$credentials$oauth_token,
    oauth_token_secret = token$credentials$oauth_token_secret
  ))
  
  if (!is.null(bodyParams)) {
    oauth <- c(oauth, bodyParams)
  }
  
  # Collect params, oauth.encode, sort and concatenated into a single string
  params <- c(urlParsed$query, oauth)
  params_esc <- setNames(oauth.encode(params), oauth.encode(names(params)))
  params_srt <- sort.names(params_esc)
  params_str <- paste0(names(params_srt), "=", params_srt, collapse = "&")
  # Generate hmac signature
  key <- paste0(oauth.encode(consumer_secret), "&", oauth.encode(token$credentials$oauth_token_secret))
  base_string <- paste0(connection, "&", oauth.encode(url), "&",
                        oauth.encode(params_str))
  oauth$oauth_signature <- hmac_sha1(key, base_string)
  sort.names(oauth)
  
  opts <- list(
    verbose = TRUE, 
    httpheader = c("Accept"="application/json", 
                   "Authorization"=paste("OAuth",
                                         paste(paste0("oauth_consumer_key=\"",oauth$oauth_consumer_key,"\""),
                                               paste0("oauth_nonce=\"",oauth$oauth_nonce,"\""),
                                               paste0("oauth_signature=\"",oauth$oauth_signature,"\""),
                                               "oauth_signature_method=\"HMAC-SHA1\"",
                                               paste0("oauth_timestamp=\"",oauth$oauth_timestamp,"\""),
                                               paste0("oauth_token=\"",oauth$oauth_token,"\""),
                                               paste0(oauth.encode(names(bodyParams)),"=\"", oauth.encode(as.character(bodyParams)), "\"", collapse=","), sep=", "),
                                         sep=" ")),
    useragent = "RCurl"
  )
  
  if(connection=="POST")
    return(postForm(url, .params=params, .opts=opts, style="post", .encoding="utf-8", .contentEncodeFun=curlPercentEncode))
  
  return(getForm(url, .params=params, .opts=opts, .encoding="utf-8", .contentEncodeFun=curlPercentEncode))
  
}
