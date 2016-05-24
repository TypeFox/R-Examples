sets_options <-
local({
    options <- list(quote = TRUE, hash = TRUE, openbounds = "()")
    function(option, value) {
        if (missing(option)) return(options)
        if (missing(value))
            options[[option]]
        else
            options[[option]] <<- value
    }
})
