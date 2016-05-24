initScore <-
function (y, offset, parms, wt) 
{
    if (!is.null(offset)) 
        y <- y - offset
    list(y = y, parms = 0, numresp = 1, numy = 1, summary = function(yval, 
        dev, wt, ylevel, digits) {
       paste("  predicted score=", format(yval, justify = "left"), 
            "  expected loss=", format(dev/wt, digits))
    }, text = function(yval, dev, wt, ylevel, digits, n, use.n) {
       
        if (use.n) {
            out <- paste(format(yval, justify = "centre"), "\n", 
                "N = ", wt, sep = "")
        } else {
            out <- format(yval, justify = "centre")
        }
        return(out)
    })
}

