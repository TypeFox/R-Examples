person_latest <-
function(api_key, page=1){
    
    url <- fromJSON(GET(url=paste("http://api.themoviedb.org/3/person/latest?api_key=", 
                                  api_key, sep=""))$url)
    
    return(url)
    
}
