#' Truncate a data frame with ellipses.
#' 
#' Prints the specified number of rows of a data frame, followed by a row of ellipses. Useful for piping to \code{knitr::kable()} for printing a truncated table in a markdown document.
#' 
#' @author Stephen Turner
#' 
#' @param df A data.frame.
#' @param n The number of rows to show before an ellipses row.
#' @import dplyr
#' 
#' @return A data frame truncated by a row of ellipses.
#' 
#' @examples
#' ellipses(mtcars, 5)
#' 
#' @export
ellipses <- function(df, n=5L) {
    stopifnot("data.frame" %in% class(df))
    els <- rep("...", ncol(df)) %>% 
        matrix(nrow=1, dimnames=list(NULL, names(df))) %>% 
        data.frame(stringsAsFactors=FALSE) %>% 
        tbl_df
    df %>% 
        head(n) %>% 
        lapply(as.character) %>% 
        data.frame(stringsAsFactors=FALSE) %>% 
        tbl_df %>% 
        bind_rows(els) %>% 
        return
}
