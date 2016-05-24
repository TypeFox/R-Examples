print.summaryStats <-
function (x, ...) 
{
    new.x <- unclass(x)
    stats.in.rows <- attr(x, "stats.in.rows")
    drop0trailing <- attr(x, "drop0trailing")
    rn <- rownames(x)
    if (stats.in.rows) {
        p.names <- c("p.value", "ChiSq_p", "Fisher_p", "Exact_p")
        index <- unlist(sapply(p.names, grep, rn))
        if (any(index)) {
            p <- new.x[index, "Combined"]
            p.char <- format(new.x)[index, "Combined"]
            if (length(grep("e", p.char))) {
                new.x[index, "Combined"] <- 0
                new.x <- format(new.x, drop0trailing = drop0trailing)
                new.x[index, "Combined"] <- p
            }
            else {
                new.x <- format(new.x, drop0trailing = drop0trailing)
            }
        }
        else {
            new.x <- format(new.x, drop0trailing = drop0trailing)
        }
        new.x[grep("NA", new.x)] <- ""
    }
    else {
        new.x <- data.frame(new.x, check.names = FALSE)
        new.x[is.na(new.x)] <- ""
    }
    print(new.x, quote = FALSE, ...)
    invisible(x)
}
