search_multi <-
function(api_key, query, page=1, include_adult=NA, language=NA){
    
    if(page<1 || page>1000){
        stop("page must be a number between 1 and 1000")
    }
    
    if(!is.logical(include_adult)){
        stop("include_adult must have a logical value.")
    }
    
    if(!is.na(language) && !is.character(language)){
        stop("language must be a ISO639-1 code")
    }
    
    l <- list(page=page, include_adult=include_adult, language=language)
    l <- l[!is.na(l)]
    
    params <- paste("&", names(l), "=",l, sep="", collapse="")
    
    url <- fromJSON(GET(URLencode(url<-paste("http://api.themoviedb.org/3/search/multi?api_key=", 
                                             api_key, "&query=", query, params, sep="")))$url)
    
    return(url)
    
}
