list_get <-
function(api_key, id){
    
    url <- fromJSON(GET(url=paste("http://api.themoviedb.org/3/list/", id, "?api_key=", 
                                  api_key, sep=""))$url)
    
    return(url)
    
}
