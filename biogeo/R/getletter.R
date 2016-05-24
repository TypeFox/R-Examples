getletter <-
function (x) 
{
    L <- rep(NA, length(x))
    S <- str_detect(x, "S")
    N <- str_detect(x, "N")
    E <- str_detect(x, "E")
    W <- str_detect(x, "W")
    L[S] <- "S"
    L[N] <- "N"
    L[E] <- "E"
    L[W] <- "W"
    x <- str_replace_all(x, "S", "")
    x <- str_replace_all(x, "N", "")
    x <- str_replace_all(x, "E", "")
    x <- str_replace_all(x, "W", "")
    x <- str_trim(x, side = "both")
    datx <- data.frame(x, L, stringsAsFactors = F)
    return(datx)
}
