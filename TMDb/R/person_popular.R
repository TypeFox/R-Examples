person_popular <-
function(api_key, page=1){

    url <- fromJSON(GET(url=paste("http://api.themoviedb.org/3/person/popular?api_key=", 
                                      api_key, "&page=", page, sep=""))$url)
    
    return(url)
    
}
