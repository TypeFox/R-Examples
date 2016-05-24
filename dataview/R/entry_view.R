#' Display vectors, lists or rows of a data frames in key-value-pairs.
#'
#' Color coded according to class of contents.
#'
#' @param x List or data frame.
#' @param i Row number to show. Press down/up to browse.
#' @return Nothing.
#' @examples
#' entry.view(Sys.getenv())
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @export
entry.view <- function(x, i=1){
    try({
        suppressMessages(attach(x)) # To get tab completion
        on.exit(detach(x))
    }, silent=TRUE)
    tryCatch({
        x.names <- if(is.null(names(x))) seq_along(x) else names(x)
        nc <- max(nchar(x.names))

        browsing <- TRUE
        while(browsing){
            cat("\n")
            cat(sprintf(sprintf("%%%is: %%s", nc),
                        x.names,
                        sapply(if(!is.null(nrow(x))) x[i,] else x,
                               function(x) style.auto(x, x))),
                sep="\n")
            if(!is.null(nrow(x)) && nrow(x) > 1){
                inchar <- readline(sprintf("\nRow %i of %i. (n)ext, (p)revious, (f)irst, (l)ast, [number], [field]: ",
                                           i, nrow(x)))
                if(inchar %in% c("n", "p", "f", "l")){
                    i <- switch(inchar, n=i+1, p=i-1, f=1, l=nrow(x))
                } else if(grepl("^\\d+$", inchar) && !inchar %in% x.names){
                    i <- as.integer(inchar)
                } else if(inchar == ""){
                    browsing <- FALSE
                } else if(inchar == "Q" && !"Q" %in% x.names){
                    interrupt()
                } else {
                    inchar <- pmatch(inchar, x.names)
                    if(!is.na(inchar)){
                        print(summary(x[[inchar]]))
                        readline("Press enter to continue.")
                    } else {
                        cat("Element not found.\n")
                    }
                }
                browsing <- browsing && 1 <= i && i <= nrow(x)
            } else {
                browsing <- FALSE
            }
        }
    }, interrupt = function(cond){
        cat(style.clear())
        signalCondition(cond)
    })
}
