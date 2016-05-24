search_movie <-
function(api_key, query, page=1, include_adult=NA, language=NA, year=NA,
                         primary_release_year=NA, search_type="phrase"){
    
    if(page<1 || page>1000){
        stop("page must be a number between 1 and 1000")
    }
    
    if(!is.logical(include_adult)){
        stop("include_adult must have a logical value.")
    }
    
    if(!is.na(language) && !is.character(language)){
        stop("language must be a ISO639-1 code")
    }
    
    if(!is.na(year) && !is.character(year)){
        stop("year must be a date string like YYYY")
    }
    
    if(!is.na(primary_release_year) && !is.character(primary_release_year)){
        stop("primary_release_year must be a date string like YYYY")
    }
    
    if(search_type!="phrase" && search_type!="ngram"){
        stop("search_type can assume only the following values: phrase, ngram")
    }
    
    l <- list(page=page, include_adult=include_adult, language=language, year=year,
              primary_release_year=primary_release_year, search_type=search_type)
    l <- l[!is.na(l)]
    
    params <- paste("&", names(l), "=",l, sep="", collapse="")
    
    url <- fromJSON(GET(URLencode(url<-paste("http://api.themoviedb.org/3/search/movie?api_key=", 
                                             api_key, "&query=", query, params, sep="")))$url)
    
    return(url)
    
}
