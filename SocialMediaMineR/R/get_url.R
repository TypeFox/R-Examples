get_url <-
function(links, sleep.time=0, return.df=F){
    final_url <- function(x) try(httr::HEAD(x)$url)
    final_url2 <- function(x) {
        Sys.sleep(sleep.time)
        cat(".")
        final_url(x)
    }
    if(return.df==F) { 
        return(vapply(links, final_url2, character(1)))
    } else {
        temp <- vapply(links, final_url2, character(1))
        url.result <- data.frame(originalURL=names(temp), resolvedURL=as.character(temp))
        return(url.result)
    }
}
