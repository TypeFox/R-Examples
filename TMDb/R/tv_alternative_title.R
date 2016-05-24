tv_alternative_title <-
function(api_key, id){

    url <- fromJSON(GET(url=paste("http://api.themoviedb.org/3/tv/", id, "/alternative_titles?api_key=", 
                                      api_key, sep=""))$url)
    
    return(url)
    
}
