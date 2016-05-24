tv_season_images <-
function(api_key, id, season_number, language=NA, include_image_language=NA){
    
    if(!is.na(language) && !is.character(language)){
        stop("language must be a ISO639-1 code")
    }
    
    if((!is.na(include_image_language) && !is.character(include_image_language)) 
       || (length(include_image_language)>5)){
        stop("include_image_language must be comma separated string with max 5 elements ISO 69-1        code.")
    } 
    
    l <- list(language=language, include_image_language=include_image_language)
    l <- l[!is.na(l)]
    
    if(length(l)>0){
        params <- paste("&", names(l), "=",l, sep="", collapse="")
        url <- fromJSON(GET(url=paste("http://api.themoviedb.org/3/tv/", id, "/season/", season_number, 
                                      "/images?api_key=", api_key, params, sep=""))$url)
    } else{
        url <- fromJSON(GET(url=paste("http://api.themoviedb.org/3/tv/", id, "/season/", season_number, 
                                      "/images?api_key=", api_key, sep=""))$url)
    }
    
    return(url)
    
}
