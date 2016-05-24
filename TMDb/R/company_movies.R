company_movies <-
function(api_key, id, page=1, language=NA, append_to_response=NA){
    
    if(page<1 || page>1000){
        stop("page must be a number between 1 and 1000")
    }
    
    if(!is.na(language) && !is.character(language)){
        stop("language must be a ISO639-1 code")
    }
    
    if(!is.na(append_to_response) && !(append_to_response %in% "movies")){
        stop("append_to_response can be NA or movies string")
    }
    
    l <- list(page=page, language=language, append_to_response=append_to_response)
    l <- l[!is.na(l)]
    
    if(length(l)>0){
        params <- paste("&", names(l), "=",l, sep="", collapse="")
        url <- fromJSON(GET(url=paste("http://api.themoviedb.org/3/company/", id, "/movies?api_key=", 
                                      api_key, params, sep=""))$url)   
    } else{
        url <- fromJSON(GET(url=paste("http://api.themoviedb.org/3/company/", id, "/movies?api_key=", 
                                      api_key, sep=""))$url)        
    }
    
    return(url)
    
}
