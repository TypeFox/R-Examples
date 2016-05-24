string.break.line <-
function (s) 
{
    if (!is.character(s)) 
        s <- as.character(s)
    string.break.line1 <- function(s) {
        chars <- substring(s, seq(len = nchar(s)), seq(len = nchar(s)))
        breaks <- c(TRUE, chars == "\n", TRUE)
        breaks <- seq(along = breaks)[breaks] - 1
        substring(s, breaks[-length(breaks)] + 1, breaks[-1] - 
            1)
    }
    lapply(s, string.break.line1)
}
