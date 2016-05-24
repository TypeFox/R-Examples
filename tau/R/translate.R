translate_Unicode_latin_ligatures <-
function(x)
{
    ## Note that we cannot use chartr() ...
    ## Note that using \uxxxx is not portable (grr) for version of R
    ## prior to 2.10.0.
    old <- Unicode_latin_ligatures_db[, 1L]
    new <- Unicode_latin_ligatures_db[, 2L]
    for(i in seq_along(old))
        x <- gsub(old[i], new[i], x, fixed = TRUE)
    x
}
