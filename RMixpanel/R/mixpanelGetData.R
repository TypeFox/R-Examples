mixpanelGetData <- function(
  account,                  # Mixpanel account
  method="export/",         # Route on Mixpanel API.
  args,                     # Argument list, e.g. list(from_date=.., to_date=.., ...)
  validitySeconds=60,
  verbose=TRUE,
  loadFromFile=FALSE,       # Load persisted file instead of download from API?
  fileName="",              # File name to persist data (defaults to './temp.txt'). If supplied, account$dataPath is used as path.
  data=FALSE,               # If FALSE, data is downloaded, but not returned as function argument.
  retryCount=0              # 
) {
  if (fileName == "") {
    filePath = paste("temp_", uuid::UUIDgenerate(), ".txt", sep="")
    on.exit( { unlink(filePath) } )
    
  } else {
    dir.create(account$dataPath, recursive=TRUE, showWarnings=FALSE)
    if(method == "export/") 
      filePath = paste(account$dataPath, "/", fileName, args$from_date, ".txt", sep="") # Events of 1 day.
    else
      filePath = paste(account$dataPath, "/", fileName, ".txt", sep="")
  }
  
  if (!loadFromFile) {
    secondsSince1970 = trunc(unclass(Sys.time()))
    endpoint = if(method == "export/") 'http://data.mixpanel.com/api/2.0' else 'http://mixpanel.com/api/2.0'
    
    args$expire = secondsSince1970 + validitySeconds
    args$api_key = account$apiKey
    args = args[order(names(args))]
    
    paramsJoined = ""
    for (name in names(args)) {
      paramsJoined = paste(paramsJoined, name, "=", args[name], sep="")
    }
    paramsJoined = paste(paramsJoined, account$apiSecret, sep='')
    hash = digest::digest(paramsJoined, "md5", serialize=FALSE)
    
    url = paste(endpoint, "/", method, "?", sep="")
    for (name in names(args)) {
      arg = args[[name]]
      if (is.character(arg))
        arg = URLencode(arg, reserved=TRUE)
      url = paste(url, name, "=", arg, "&", sep="")
    }
    url = paste(url, "sig=", hash, sep="")
    if (verbose)
      cat("## Download ", url, "...\n", sep="")
    
    for (trial in 0:retryCount) {
      unlink(filePath)
      code = download.file(url, filePath, "curl", quiet=TRUE)
      if (code == 0) break
      cat("## !!! Retry... ", date(), "\n")
      Sys.sleep(2)
    }
  }
  
  if (data) {
    if (verbose)
      cat("## Scan... ", date(), "\n")
    res = scan(filePath, what=character(), sep="\n", quiet=!verbose)
    res
  }
}
