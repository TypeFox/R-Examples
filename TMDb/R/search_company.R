search_company <-
function(api_key, query, page=1){
    
    if(page<1 || page>1000){
        stop("page must be a number between 1 and 1000")
    }
    
    url <- fromJSON(GET(URLencode(url<-paste("http://api.themoviedb.org/3/search/company?api_key=", 
                                  api_key, "&query=", query, "&page=", page, sep="")))$url)
    
    return(url)
    
}
