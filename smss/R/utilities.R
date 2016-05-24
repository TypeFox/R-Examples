describe_df <- function(x) {
  sprintf("\\code{data.frame} with %d observations of the following %d variables,",
          nrow(x), ncol(x))
}

df_format <- function(x) {
  template <- "%s\n\\describe{\n%s\n}"
  items <- sapply(names(x),
                  function(i) sprintf("\\item{\\code{%s}}{%s}", i,
                                      ifelse(is.null(comment(x[[i]])), i, comment(x[[i]]))))
  sprintf(template, describe_df(x), paste(items, collapse="\n"))
}
