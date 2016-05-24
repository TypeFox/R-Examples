tv_airing_today <-
function(api_key, page=1, language=NA, timezone=NA){
    
    if(page<1 || page>1000){
        stop("page must be a number between 1 and 1000")
    }
    
    if(!is.na(language) && !is.character(language)){
        stop("language must be a ISO639-1 code")
    }
    
    l <- list(page=page, language=language, timezone=timezone)
    l <- l[!is.na(l)]
    
    if(length(l)>0){
        params <- paste("&", names(l), "=",l, sep="", collapse="")
        url <- fromJSON(GET(url=paste("http://api.themoviedb.org/3/tv/airing_today?api_key=", 
                                      api_key, params, sep=""))$url)
    } else{
        url <- fromJSON(GET(url=paste("http://api.themoviedb.org/3/tv/airing_today?api_key=", 
                                      api_key, sep=""))$url)
    }
    
    return(url)
    
}
