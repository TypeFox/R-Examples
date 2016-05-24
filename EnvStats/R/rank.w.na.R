rank.w.na <-
function (x, na.action = c("remove", "mid.rank", "lowest", "highest", 
    "include"), warn = TRUE) 
{
    if ((n.na <- sum(na.index <- is.na(x))) == 0) 
        ranks <- rank(x)
    else {
        na.action <- match.arg(na.action)
        if (na.action == "include") {
            ranks <- rank(x)
            ranks[is.na(x)] <- NA
            if (warn) {
                string1 <- ifelse(n.na > 1, "values", "value")
                warning(paste(n.na, "missing", string1, "ignored in ranking"))
            }
        }
        else if (na.action != "mid.rank") {
            na.last <- switch(na.action, lowest = FALSE, highest = TRUE, 
                remove = NA)
            ranks <- rank(x, na.last = na.last)
            if (warn) {
                string1 <- ifelse(n.na > 1, "values", "value")
                string2 <- ifelse(n.na > 1, "ranks", "rank")
                warn <- switch(na.action, lowest = paste(n.na, 
                  "missing", string1, "assigned lowest", string2), 
                  highest = paste(n.na, "missing", string1, "assigned highest", 
                    string2), remove = paste(n.na, "missing", 
                    string1, "removed before ranking"))
                warning(warn)
            }
        }
        else {
            ranks <- numeric(length(x))
            ranks[!na.index] <- rank(x[!na.index])
            ranks[na.index] <- mean(ranks[!na.index])
            if (warn) {
                string <- ifelse(n.na > 1, "values", "value")
                warning(paste(n.na, "missing", string, "set to the mean rank of the", 
                  "non-missing values"))
            }
        }
    }
    ranks
}
