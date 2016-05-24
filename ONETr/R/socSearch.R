socSearch <-
function(soc){
    # call api and parse xml
    baseURL <- "https://services.onetcenter.org/ws/online/occupations/"
    output <- getURL(paste(baseURL,soc,"/details",sep=""), userpwd=paste(get("creds",envir=cacheEnv)[[1]],":",get("creds",envir=cacheEnv)[[2]], sep=""), httpauth = 1L)
    if(grepl("Authorization Error",output)){
       message("Error: Your API credentials are invalid. Please enter valid HTTPS credentials using setCreds().")
    }
    else if(grepl("<error>",output)){
       message("Error: A valid O*NET-SOC code is required.")
    }
    else{
      raw <- xmlInternalTreeParse(output)
      rootNode <- xmlRoot(raw)
      list <- xmlToList(rootNode)
      return(list)
    }
}
