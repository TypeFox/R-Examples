# Copied from seewave?

ftwindow <-
function (wl, wn = "hanning") 
{
    if (wn == "bartlett") 
        w <- bartlett.w(wl)
    if (wn == "blackman") 
        w <- blackman.w(wl)
    if (wn == "flattop") 
        w <- flattop.w(wl)
    if (wn == "hamming") 
        w <- hamming.w(wl)
    if (wn == "hanning") 
        w <- hanning.w(wl)
    if (wn == "rectangle") 
        w <- rectangle.w(wl)
    return(w)
}
