#' Search data dictionary.
#'
#' This function is used to search the College Scorecard
#' data dictionary.
#'
#' @param search_string Character string for search. Can use regular
#'     expression for search. Must escape special characters,
#'     \code{. \ | ( ) [ \{ ^ $ * + ?}, with a doublebackslash
#'     \code{\\\\}.
#' @param search_col Column to search. The default column is the
#'     "description" column. Other options include: "varname",
#'     "dev_friendly_name", "dev_category", "label".
#' @param ignore_case Search is case insensitive by default. Change to
#'     \code{FALSE} to restrict search to exact case matches.
#' @param limit Only the first 10 dictionary items are returned by
#'     default. Increase to return more values. Set to \code{Inf} to
#'     return all items matched in search'
#' @param confirm Use to confirm status of variable name in
#'     dictionary. Returns \code{TRUE} or \code{FALSE}.
#'
#' @examples
#' sc_dict('state')
#' sc_dict('^st', search_col = 'varname') # variable names starting with 'st'
#' sc_dict('.', limit = Inf) # return full dictionary (not recommended)

#' @export
sc_dict <- function(search_string,
                    search_col = c('description',
                                   'varname',
                                   'dev_friendly_name',
                                   'dev_category',
                                   'label'),
                    ignore_case = TRUE,
                    limit = 10,
                    confirm = FALSE) {

    ## only for confirm
    if (confirm) {
        return(!is.null(hash[[search_string]]))
    }

    ## get values
    rows <- grepl(search_string,
                  dict[[match.arg(search_col)]],
                  ignore.case = ignore_case)

    if (all(rows == FALSE)) {
        return(cat('\nNo matches! Try again with new string or column.\n\n'))
    }

    ## pull data
    out <- dict[rows,]

    ## get unique varnames
    uniqv <- unique(out[['varname']])

    ## pretty print
    for (i in 1:min(length(uniqv), limit)) {

        ## subset
        d <- out[out[['varname']] == uniqv[i],]

        ## console table
        cat('\n' %+% paste(rep('', 70), collapse = '-') %+% '\n')
        cat('varname: ' %+% d[['varname']][1])
        cat(rep('', 51 - nchar(d[['varname']][1]) -
                    nchar(d[['dev_category']][1])))
        cat('category: ' %+% d[['dev_category']][1])
        cat('\n' %+% paste(rep('', 70), collapse = '-') %+% '\n')
        cat('DEVELOPER FRIENDLY NAME:\n\n')
        cat(d[['developer_friendly_name']][1] %+% '\n\n')
        cat('DESCRIPTION:\n\n')
        cat(strwrap(d[['description']][1], 70) %+% '\n')
        cat('\n')
        cat('VALUES: ')
        if (is.na(d[['value']][1])) {
            cat('NA\n\n')
        } else {
            cat('\n\n')
            for (j in seq(nrow(d))) {

                cat(d[['value']][j] %+% ' = ' %+% d[['label']][j] %+% '\n')

            }
            cat('\n')
        }

    }

    cat('Printed information for ' %+% min(length(uniqv), limit) %+% ' of out ')
    cat(length(uniqv) %+% ' variables.\n')
    if (limit < length(uniqv)) cat('Increase limit to see more variables.\n')
    cat('\n')
}
