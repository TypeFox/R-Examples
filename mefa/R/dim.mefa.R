`dim.mefa` <-
function (x) 
{
    nseg <- if (is.null(x$segm)) 1 else length(x$segm)
    out <- c(nrow(x$xtab), ncol(x$xtab), nseg)
    return(out)
}

