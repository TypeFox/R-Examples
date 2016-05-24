movie_popular <-
function(api_key, page=1, language=NA){
    
    l <- list(page=page, language=language)
    l <- l[!is.na(l)]
    
    if(length(l)>0){
        params <- paste("&", names(l), "=",l, sep="", collapse="")
        url <- fromJSON(GET(url=paste("http://api.themoviedb.org/3/movie/popular?api_key=", 
                                      api_key, params, sep=""))$url)
    } else{
        url <- fromJSON(GET(url=paste("http://api.themoviedb.org/3/movie/popular?api_key=", 
                                      api_key, sep=""))$url)
    }
    
    return(url)
    
}
