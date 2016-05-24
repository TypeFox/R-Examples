tv_content_ratings <-
function(api_key, id){
    
    url <- fromJSON(GET(url=paste("http://api.themoviedb.org/3/tv/", id, "/content_ratings?api_key=", 
                                  api_key, sep=""))$url)
    
    return(url)
    
}
