search_keyword <-
function(api_key, query, page=1){
    
    if(page<1 || page>1000){
        stop("page must be a number between 1 and 1000")
    }
    
    l <- list(page=page)
    l <- l[!is.na(l)]
    
    params <- paste("&", names(l), "=",l, sep="", collapse="")
    url <- fromJSON(GET(URLencode(url<-paste("http://api.themoviedb.org/3/search/collection?api_key=", 
                                             api_key, "&query=", query, params, sep="")))$url)
    
    return(url)
    
}
