stop.redef <-
function (locstring = "", ...)
{
    print(locstring, quote = FALSE)
    base::stop(locstring, ...)
}
