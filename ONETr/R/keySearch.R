keySearch <-
function(keyword){
      output <- getURL(paste("https://services.onetcenter.org/ws/mnm/search?keyword=",keyword,sep=""), userpwd=paste(get("creds",envir=cacheEnv)[[1]],":",get("creds",envir=cacheEnv)[[2]], sep=""), httpauth = 1L)
      if(grepl("Authorization Error",output)){
        message("Your API credentials are invalid. Please enter valid HTTPS credentials using setCreds().")
      }
      else if(grepl('total="0"',output)){
        message("Your keyword returned no results. Please try another keyword or occupational title.")
      }
      else{
        output <- xmlParse(output)
        keyOutput <- xmlToDataFrame(nodes = getNodeSet(output, "//career"))
        message("Find a SOC code below and use with socSearch() function to pull job data.")
        return(keyOutput[,1:2])
      }
      
}
