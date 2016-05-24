number.suffix.scalar <-
function (x) 
{
    if (!(length(x) == 1) || !is.numeric(x) || is.na(x)) 
        stop("'x' must be a numeric scalar")
    x <- as.character(x)
    n <- nchar(x)
    last.digit <- substring(x, n, n)
    switch(last.digit, `0` = "'th", `1` = "'st", `2` = "'nd", 
        `3` = "'rd", `4` = "'th", `5` = "'th", `6` = "'th", `7` = "'th", 
        `8` = "'th", `9` = "'th")
}
