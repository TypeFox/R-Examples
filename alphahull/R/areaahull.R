areaahull <-
function (x, timeout = 5) 
{
    area <- evalWithTimeout(try(areaahulleval(x), silent = TRUE), 
        timeout = timeout)
    if (!is.numeric(area)) {
        warning("Problem in area computation (Returns NA)")
        area <- NA
    }
    if (is.numeric(area) & area < 0) {
        warning("Problem in area computation (Returns NA)")
        area <- NA
    }
    return(area)
}
