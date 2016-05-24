person_external_ids <-
function(api_key, id){
    
    url <- fromJSON(GET(url=paste("http://api.themoviedb.org/3/person/", id, "/external_ids?api_key=", 
                                  api_key, sep=""))$url)
    
    return(url)
    
}
