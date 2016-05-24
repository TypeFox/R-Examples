search_tv <-
function(api_key, query, page=1, language=NA, 
                      first_air_date_year=NA, search_type="phrase"){
    
    if(page<1 || page>1000){
        stop("page must be a number between 1 and 1000")
    }
    
    if(!is.na(language) && !is.character(language)){
        stop("language must be a ISO639-1 code")
    }
    
    if(!is.na(first_air_date_year) && !is.character(first_air_date_year)){
        stop("first_air_date_year must be a date string like YYYY")
    }
    
    if(search_type!="phrase" && search_type!="ngram"){
        stop("search_type can assume only the following values: phrase, ngram")
    }
    
    l <- list(page=page, first_air_date_year=first_air_date_year, search_type=search_type)
    l <- l[!is.na(l)]
    
    params <- paste("&", names(l), "=",l, sep="", collapse="")
    
    url <- fromJSON(GET(URLencode(url<-paste("http://api.themoviedb.org/3/search/tv?api_key=", 
                                             api_key, "&query=", query, params, sep="")))$url)
    
    return(url)
    
}
