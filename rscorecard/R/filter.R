#' Filter scorecard data by variable values.
#'
#' This function is used to filter the downloaded scorecard data. It
#' converts idiomatic R into the format required by the API call.
#'
#' @param sccall Current list of parameters carried forward from prior
#'     functions in the chain (ignore)
#' @param ... Expressions to evaluate
#'
#' @examples
#' \dontrun{
#' sc_filter(region == 1) # New England institutions
#' sc_filter(stabbr == c('TN','KY')) # institutions in Tennessee and Kentucky
#' sc_filter(control != 3) # exclude private, for-profit institutions
#' sc_filter(control == c(1,2)) # same as above
#' sc_filter(control == 1:2) # same as above
#' sc_filter(stabbr == 'TN', control == 1, locale == 41:43) # TN rural publics
#' }

#' @export
sc_filter <- function(sccall, ...) {

    ## check first argument
    if (identical(class(try(sccall, silent = TRUE)), 'try-error')) {
        stop('Chain not properly initialized. Be sure to start with sc_init().',
             call. = FALSE)
    }

    ## get expressions
    filter <- lapply(lazyeval::lazy_dots(...),
                     function(x) as.list(bquote(.(x[['expr']]))))

    ## error handling
    for (i in 1:length(filter)) {
        if (!identical(filter[[i]][[1]], as.symbol('=='))
            && !identical(filter[[i]][[1]], as.symbol('!='))) {
            stop('Must use either \"==\" or \"!=\" in sc_filter.',
                 call. = FALSE)
        }
        if (!sc_dict(tolower(as.character(filter[[i]][[2]])), confirm = TRUE)) {
            stop('Variable \"' %+% filter[[i]][[2]]
                 %+% '\" not found in dictionary. '
                 %+% 'Please check your spelling or search dictionary: '
                 %+% '?sc_dict()', call. = FALSE)
        }
    }

    ## convert to developer-friendly names
    if (!sccall[['dfvars']]) {
        for (i in 1:length(filter)) {
            filter[[i]][[2]] <- hash[[tolower(as.character(filter[[i]][[2]]))]]
        }
    }

    ## grab categories
    cats <- vapply(filter, function(x) { hash[[as.character(x[[2]]) %+% '_c']] },
                   character(1), USE.NAMES = FALSE)

    ## convert idiomatic R to scorecard API style
    filter <- vapply(filter, api_convert, character(1))

    ## paste, clean, and return
    filter <- paste(cats %+% '.' %+% filter, collapse = '&')
    filter <- gsub('root.', '', filter, fixed = TRUE)
    sccall[['filter']] <- filter
    sccall

}



