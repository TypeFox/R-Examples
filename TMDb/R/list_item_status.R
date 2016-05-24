list_item_status <-
function(api_key, id, movie_id){
    
    url <- fromJSON(GET(url=paste("http://api.themoviedb.org/3/list/", id, "/item_status?api_key=", 
                                  api_key, "&movie_id=", movie_id, sep=""))$url)
    
    return(url)
    
}
