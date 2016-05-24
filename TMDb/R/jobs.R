jobs <-
function(api_key){

    url <- fromJSON(GET(url=paste("http://api.themoviedb.org/3/job/list?api_key=", 
                                      api_key, sep=""))$url)
    
    return(url)
    
}
