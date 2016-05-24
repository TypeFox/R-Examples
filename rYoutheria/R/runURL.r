#' @import plyr
 
runURL <-
  function(URL='', type=NULL){
    
    if(is.null(URL)) stop('URL is NULL')
    if(is.null(type)) stop('type is NULL')
    
    apiURL <- 'http://www.utheria.org/api/'
    
    if(type=='m') baseURL <- paste(apiURL, 'ValueByType', sep = '')
    if(type=='l') baseURL <- paste(apiURL, 'Location', sep = '')
    if(type=='t') baseURL <- paste(apiURL, 'Measurement', sep = '')
    if(type=='c') baseURL <- paste(apiURL, 'Country', sep = '')
    
    fullURL <- paste(baseURL,URL,sep='')
    # replace spaces in the URL with the html '%20'
    fullURL <- gsub(' ', '%20', fullURL)
    
    if(length(fullURL) == 1){
      
      out <- executeURL(fullURL, type)
    
    } else {
      
      outMulti <- lapply(fullURL, executeURL, type)
      out <- ldply(outMulti, data.frame, stringsAsFactors=FALSE)
      
    }
        
    return(out)
    
  }