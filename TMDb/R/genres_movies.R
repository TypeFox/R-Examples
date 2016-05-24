genres_movies <-
function(api_key, id, page=1, language=NA, include_all_movies=NA, include_adult=NA){
    
    if(page<1 || page>1000){
        stop("page must be a number between 1 and 1000")
    }
    
    if(!is.na(language) && !is.character(language)){
        stop("language must be a ISO639-1 code")
    }
    
    if(!is.na(include_all_movies) && !is.logical(include_all_movies)){
        stop("include_all_movies must have a logical value.")
    }
    
    if(!is.na(include_adult) && !is.logical(include_adult)){
        stop("include_adult must have a logical value.")
    }
    
    l <- list(page=page, language=language, include_all_movies=include_all_movies, include_adult=include_adult)
    l <- l[!is.na(l)]
    
    if(length(l)>0){
        params <- paste("&", names(l), "=",l, sep="", collapse="")
        url <- fromJSON(GET(url=paste("http://api.themoviedb.org/3/genre/", id, "/movies?api_key=", 
                                      api_key, params, sep=""))$url)   
    } else{
        url <- fromJSON(GET(url=paste("http://api.themoviedb.org/3/genre/", id, "/movies?api_key=", 
                                      api_key, sep=""))$url)        
    }
    
    return(url)
    
}
