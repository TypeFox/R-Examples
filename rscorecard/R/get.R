#' Get scorecard data.
#'
#' This function gets the College Scorecard data by compiling and
#' converting all the previous piped output into a single URL string
#' that is used to get the data.
#'
#' @param sccall Current list of parameters carried forward from prior
#'     functions in the chain (ignore)
#' @param api_key Personal API key requested from
#'     \url{https://api.data.gov/signup} stored in a string. If you
#'     first set your key using \code{sc_key}, then you may omit this
#'     parameter. A key set here will take precedence over any set in
#'     the environment (DATAGOV_API_KEY).
#'
#' @examples
#' \dontrun{
#' sc_get('<API KEY IN STRING>')
#' key <- '<API KEY IN STRING>'
#' sc_get(key)
#' }
#'
#' @section Obtain a key:
#' To obtain an API key, visit \url{https://api.data.gov/signup}

#' @export
sc_get <- function(sccall, api_key) {

    ## check first argument
    if (identical(class(try(sccall, silent = TRUE)), 'try-error')) {
        stop('Chain not properly initialized. Be sure to start with sc_init().',
             call. = FALSE)
    }

    ## add year
    re <- '(=|,)(' %+% 'academics' %+|%
                'admissions' %+|%
                'aid' %+|%
                'completion' %+|%
                'cost' %+|%
                'earnings' %+|%
                'repayment' %+|%
                'student' %+% '\\.)'
    sccall[['select']] <- gsub(re, '\\1' %+% sccall[['year']] %+%
                                      '.' %+% '\\2', sccall[['select']])
    sccall[['filter']] <- gsub(re, '\\1' %+% sccall[['year']] %+%
                                      '.' %+% '\\2', sccall[['filter']])

    ## create url for call
    url <- 'https://api.data.gov/ed/collegescorecard/v1/schools.json?' %+%
         sccall[['filter']] %+% sccall[['select']] %+% sccall[['zip']]

    ## check for key
    if (missing(api_key)) {
        api_key <- Sys.getenv('DATAGOV_API_KEY')
        if (identical(api_key, '')) {
            stop('Missing API key; ?sc_key for details', call. = FALSE)
        }
    }

    ## first GET
    con <- url %+% '&_page=0&_per_page=100&api_key=' %+% api_key
    init <- fromJSON(con)

    ## return if no options
    if (init[['metadata']][['total']] == 0) {
        stop('No results! Broaden your search or try different variables.',
             call. = FALSE)
    }

    if (init[['metadata']][['total']] > nrow(init[['results']])) {

        ## get number of pages needed
        pages <- ceiling(init[['metadata']][['total']] / 100)

        message('Large request will require: ' %+% pages %+% ' additional pulls.')

        ## download data in chunks and bind
        page_list <- vector('list', pages)
        for (i in 1:pages) {
            message('Request chunk ' %+% i)
            con <- url %+% '&_page=' %+% i %+% '&_per_page=100&api_key=' %+% api_key
            page_list[[i]] <- fromJSON(con)[['results']]
        }

        df <- dplyr::rbind_all(page_list)
        df <- dplyr::rbind_list(init[['results']], df)

    } else {

        df <- dplyr::tbl_df(init[['results']])

    }

    ## convert names back to non-developer-friendly names and return
    if (!sccall[['dfvars']]) {
        names(df) <- vapply(names(df), function(x) {
            re <- '^[0-9]{0,4}\\.?(' %+% 'academics' %+|%
                         'admissions' %+|%
                         'aid' %+|%
                         'completion' %+|%
                         'cost' %+|%
                         'earnings' %+|%
                         'repayment' %+|%
                         'root' %+|%
                         'school' %+|%
                         'student' %+% ')\\.'
            hash[[gsub(re, '', x)]]},
            character(1), USE.NAMES = FALSE
            )
    }

    ## add year column and return
    df[['year']] <- sccall[['year']]
    df

}



