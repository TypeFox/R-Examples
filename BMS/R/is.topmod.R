is.topmod <-
function (tmo) 
{
    if (is.element("topmod", class(tmo))) 
        return(TRUE)
    else return(FALSE)
}
