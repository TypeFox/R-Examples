
api_convert <- function(filters) {
    ## convert rhs
    filters[3] <- gsub('c\\((.+)\\)', '\\1', filters[3])
    filters[3] <- gsub('\\"|\\\'', '', filters[3])
    filters[3] <- gsub(', ', ',', filters[3])
    filters[3] <- gsub(' ', '%20', filters[3])
    filters[3] <- gsub(':', '..', filters[3])
    ## convert operator
    filters[1] <- as.character(filters[[1]])
    filters[1] <- gsub('==', '=', filters[1])
    filters[1] <- gsub('!=', '__not=', filters[1])
    filters[1] <- ifelse(grepl('..', filters[3], fixed = TRUE),
                        gsub('=', '__range=', filters[1]),
                        filters[1])
    ## paste
    filters[[2]] %+% filters[[1]] %+% filters[[3]]
}
