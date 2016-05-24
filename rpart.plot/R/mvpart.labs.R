# mvpart.labs.R: generate labels for mvpart (mrt) objects

# This extends the text() function in rpart.mvpart.
# Note that extra=0 and extra=1 are compatible with that function.
# The following table shows allowable values for extra (and as always,
# add 100 to include the percent).  My personal favorites are 107 and 111.
#
#    0   dev
#    1   dev, n
#    2   dev, yhat
#    3   dev, yhat / sum(yhat)
#    4   sqrt(dev)
#    5   sqrt(dev), n
#    6   sqrt(dev), yhat
#    7   sqrt(dev), yhat / sum(yhat)
#    8   predominant species
#    9   predominant species, n
#    10  predominant species, yhat
#    11  predominant species, yhat / sum(yhat)

get.mvpart.labs <- function(x, extra, under, digits, xsep, varlen)
{
    frame <- x$frame
    yval2 <- frame$yval2 # fit per species i.e. per column of response matrix y
    nspecies <- ncol(yval2)
    if(is.null(xsep))
        xsep <- "  " # two spaces
    ex <- if(extra < 100) extra else extra - 100
    if(ex <= 3)
        main <- formate(frame$dev, digits=digits)
    else if(ex <= 7)
        main <- formate(sqrt(frame$dev), digits=digits)
    else {
        stopifnot(all(!is.na(yval2))) # needed because which.max discards NAs
        main <- apply(yval2, 1, which.max) # index of max in each row
        # convert species number to species name, if the names are available
        if(length(colnames(x$y)) == nspecies)
            main <- colnames(x$y)[main]
        main <- as.character(main)
    }
    if(ex == 3 || ex == 7 || ex == 11) # divide each row by its sum
        for (i in 1:nrow(yval2))
            yval2[i,] <- yval2[i,] / sum(yval2[i,])
    resp.per.species <- formatf(yval2, digits, strip.leading.zeros=TRUE) # TODO use formate?
    resp.per.species <- apply(matrix(resp.per.species, ncol=nspecies),
                              1, paste.with.breaks, collapse=xsep)
    newline <- if(under) "\n\n" else "\n"
    labs <-
        if(ex == 0 || ex == 4 || ex == 8)
            main
        else if(ex == 1 || ex == 5 || ex == 9)
            sprintf("%s%sn=%s", main, newline, format0(frame$n, digits))
        else if(ex == 2 || ex == 3 || ex == 6 || ex == 7 || ex == 10 || ex == 11)
            paste0(main, newline, resp.per.species)
        else
            stop0("extra=", extra, " is illegal (for method=\"", x$method, "\")")

    if(extra >= 100) { # add percent?
        sep <- switch(ex+1,   # figure out where to put percent (same line? below? etc.)
                      newline,  # 0 may be a double newline
                      "  ",     # 1
                      "\n",     # 2
                      "\n",     # 3
                      newline,  # 4
                      "  ",     # 5
                      "\n",     # 6
                      "\n",     # 7
                      newline,  # 8
                      "  ",     # 9
                      "\n",     # 10
                      "\n")     # 11

        labs <- sprintf("%s%s%s%%", labs, sep,
                        formatf(100 * frame$wt / frame$wt[1], digits=max(0, digits-2)))
    }
    labs
}
