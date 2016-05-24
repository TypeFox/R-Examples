## if you use the example configuration supplied in "code"
## then this script will be loaded for all requests
## before the run() function is evaluated. You can use it
## for global processing, for example, to handle cookies.

getCookies <- function() {
    ## raw.cookies is a variable populated by the FastRWeb engine
    ckv <- gsub("^ +","",strsplit(request$raw.cookies,";",fixed=TRUE)[[1]])
    cookies <<- if (length(ckv)) {
            keys = unlist(lapply(strsplit(ckv,"="),function(x) x[1]))
                vals = substr(ckv,nchar(keys)+2,99999)
                names(vals) = keys
                vals
          } else character(0)
}

## just a dummy example - this will create a "cookes" variable
## that can be used by all run() scripts to access cookie contents.
getCookies()
