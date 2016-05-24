cbind.no.warn <-
function (..., deparse.level = 1) 
{
    oldopts <- options(warn = -1)
    on.exit(options(oldopts))
    base::cbind(..., deparse.level = deparse.level)
}
