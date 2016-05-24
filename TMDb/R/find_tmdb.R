find_tmdb <-
function(api_key, id, external_source, language=NA){
    
    if(!is.na(language) && !is.character(language)){
        stop("language must be a ISO639-1 code")
    }
    
    l <- list(language=language)
    l <- l[!is.na(l)]
    
    if(length(l)>0){
        params <- paste("&", names(l), "=",l, sep="", collapse="")
        url <- fromJSON(GET(url=paste("http://api.themoviedb.org/3/find/", id, "?api_key=", 
                                      api_key, "&external_source=", external_source, params, sep=""))$url)   
    } else{
        url <- fromJSON(GET(url=paste("http://api.themoviedb.org/3/find/", id, "?api_key=", 
                                      api_key, "&external_source=", external_source, sep=""))$url)        
    }
    
    return(url)
    
}
