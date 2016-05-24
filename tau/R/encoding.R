
## Encoding utilities.
##
## ceeboo 2008

is.ascii <-
function(x)
    .Call(R_isASCII, x)

is.utf8 <-
function(x)
    .Call(R_isUTF8, x)

## TODO consider a fixWorkspace function

fixEncoding <-
function(x, latin1 = FALSE)
    .Call(R_fixEncoding, x, latin1)

is.locale <-
function(x)
    .Call(R_isLocale, x)

translate <-
function(x, recursive = FALSE, internal = FALSE)
{
    if (!is.null(names(x)))
        names(x) <- translate(names(x))
    if (is.character(x)) {
        if (internal) {
            x <- lapply(x, function(x) {
                if (Encoding(x) != "unknown") {
                    y <- iconv(x, from = Encoding(x))
                    if (!is.na(y))
                        x <- y
                }
                x
            })
            as.character(x)
        } else
            .Call(R_translateToLocale, x)
    }
    else if (is.factor(x))
        levels(x) <- translate(levels(x))
    else if (is.list(x) && recursive)
        lapply(x, translate, TRUE)
    else
        x
}

##
