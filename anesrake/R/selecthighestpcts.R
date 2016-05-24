selecthighestpcts <-
function(dif, inputter, pctlim, 
    tostop = 1, warn = 1) {
    if (pctlim > 1) {
        pctlim <- pctlim/100
    }
    found <- inputter[dif >= pctlim]
    if (length(found) == 0 & tostop == 1) {
        stop(paste("No variables are off by more than", pctlim * 
            100, "percent using the method you have chosen, either weighting is unnecessary or a smaller pre-raking limit should be chosen."))
    }
    if (length(found) == 0 & tostop == 0 & warn == 1) {
        warning(paste("No variables are off by more than", pctlim * 
            100, "percent using the method you have chosen, either weighting is unnecessary or a smaller pre-raking limit should be chosen."))
        print(paste("No variables are off by more than", pctlim * 
            100, "percent using the method you have chosen, either weighting is unnecessary or a smaller pre-raking limit should be chosen."))
    }
    found
}

