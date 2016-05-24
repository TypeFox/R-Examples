search_person <-
function(api_key, query, page=1, include_adult=NA, search_type="phrase"){
    
    if(page<1 || page>1000){
        stop("page must be a number between 1 and 1000")
    }
    
    if(!is.logical(include_adult)){
        stop("include_adult must have a logical value.")
    }
    if(search_type!="phrase" && search_type!="ngram"){
        stop("search_type can assume only the following values: phrase, ngram")
    }
    
    l <- list(page=page, include_adult=include_adult, search_type=search_type)
    l <- l[!is.na(l)]
    
    params <- paste("&", names(l), "=",l, sep="", collapse="")
    
    url <- fromJSON(GET(URLencode(url<-paste("http://api.themoviedb.org/3/search/person?api_key=", 
                                             api_key, "&query=", query, params, sep="")))$url)
    
    return(url)
    
}
