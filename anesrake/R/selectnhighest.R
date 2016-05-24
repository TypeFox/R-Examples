selectnhighest <-
function(dif, inputter, nlim) {
    rnk <- rank(dif)
    if (length(nlim) > length(inputter)) {
        stop(paste("cannot rake on", nlim, "variables.  Too few targets specified"))
    }
    found <- inputter[rnk <= nlim]
    found
}

