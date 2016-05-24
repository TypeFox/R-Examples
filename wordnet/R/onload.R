dict <-
local({
    d <- NULL
    function(new, ...) {
        if (!missing(new))
            d <<- new
        else
            d
    }
})

.onLoad <-
function(libname, pkgname)
{
    .jpackage(pkgname, lib.loc = libname)

    if (initDict())
        dict(getDictInstance())
}
