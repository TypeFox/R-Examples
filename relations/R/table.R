relation_table <-
function(x, memberships = TRUE)
{
    if(!is.relation(x))
        stop("Argument 'x' must be a relation.")
    DF <- as.data.frame(x)
    M <- attr(DF, "memberships")
    if (memberships && !is.null(M))
        DF <- cbind(DF, " " = paste("[", M, "]", sep = ""))
    .structure(DF, class = c("relation_table", "data.frame"))
}

print.relation_table <-
function(x, ...)
{
    y <- as.data.frame(x)
    if(!length(row.names(y)))
        print.data.frame(y)
    else {
        y <- as.matrix(format(y, justify = "left"))
        rownames(y) <- rep.int("", nrow(y))
        print(y, quote = FALSE)
    }
    invisible(x)

}
