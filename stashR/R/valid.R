setValidity("remoteDB",
            function(object) {
                if(length(object@name) == 0)
                    return("object needs a non-empty 'name'")
                if(length(object@url) == 0)
                    return("object needs a non-empty 'url'")
                if(!validURL(object@url))
                    return("object needs a 'url' of type 'http://'")
                TRUE
            })

validURL <- function(URL) {
    if(length(grep("^http://", URL)) > 0)
        TRUE
    else
        FALSE
}
