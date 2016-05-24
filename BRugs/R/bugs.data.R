"bugsData" <- 
    function(data, fileName = file.path(tempdir(), "data.txt"), format="E", digits = 5){
        if (is.character(unlist(data))) { 
            if(is.R()) {
                data.list <- lapply(as.list(data), get, pos = parent.frame(2))
                names(data.list) <- as.list(data)
                write.datafile(lapply(data.list, formatC, digits = digits, format = format), fileName)
            }
            else {
                data.list <- lapply(as.list(data), get, where = parent.frame(2))
                names(data.list) <- unlist(data)
                writeDatafileS4(data.list, towhere = "data.txt")
            }
        }
        else if(is.list(data)) { 
            data <- lapply(data, function(x){x <- if(is.character(x)||is.factor(x)) match(x, unique(x)) else x})
            if(is.R()) {
                write.datafile(lapply(data, formatC, digits = digits, format = format), fileName)
            }
            else {
                writeDatafileS4(data, towhere = "data.txt")
            }
        }
        else stop("Expected a list of data, a list or vector of variable names")      
        invisible(fileName)
    }


if(is.R()){
    ## need some fake functions for codetools
    toSingleS4 <- function(...)
      stop("This function is not intended to be called in R!")
    "writeDatafileS4" <- toSingleS4
} else {

### The rest of this file is for S-PLUS only...


"writeDatafileS4" <- 
#
# Writes to file "towhere" text defining a list containing "DATA" in a form compatable with WinBUGS.
# Required arguments:
# DATA - either a data frame or else a list consisting of any combination of scalars, vectors, arrays or data frames (but not lists).
#   If a list, all list elements that are not data.frames must be named. Names of data.frames in DATA are ignored.
# Optional arguments:
# towhere - file to receive output. Is clipboard by default, which is convenient for pasting into a WinBUTS ODC file.
# fill - If numeric, number of columns for output. If FALSE, output will be on one line. If TRUE (default), number of
#   columns is given by .Options$width.
# Value:
# Text defining a list is output to file "towhere". 
# Details:
#  The function performs considerable checking of DATA argument. Since WinBUGS requires numeric input, no factors or character vectors
# are allowed. All data must be named, either as named elements of DATA (if it is a list) or else using the names given in data frames.
# Data frames may contain matrices. 
# Arrays of any dimension are rearranged to be in row-major order, as required by WinBUGS. Scientific notation is also handled properly.
# In particular, the number will consist of a mantissa _containing a decimal point_ followed by "E", then either "+" or "-", and finally 
# a _two-digit_ number. S-Plus does not always provide a decimal point in the mantissa, uses "e" instead of "E", followed by
# either a "+" or "-" and then _three_ digits.
# Written by Terry Elrod. Disclaimer: This function is used at the user's own risk. 
# Please send comments to Terry.Elrod@UAlberta.ca.
# Revision history: 2002-11-19. Fixed to handle missing values properly.
function(DATA, towhere = "clipboard", fill = TRUE)
{
    formatDataS4 = 
    #
    # Prepared DATA for input to WinBUGS.
    function(DATA)
    {
        if(!is.list(DATA))
            stop("DATA must be a named list or data frame.")
        dlnames <- names(DATA)
        if(is.data.frame(DATA))
            DATA <- as.list(DATA)
        #
        # Checking for lists in DATA....
        lind <- sapply(DATA, is.list)
        # Checking for data frames in DATA....
        dfind <- sapply(DATA, is.data.frame)
        # Any lists that are not data frames?...
        if(any(lind & !dfind)) stop("DATA may not contain lists.")
        # Checking for unnamed elements of list that are not data frames....
        if(any(dlnames[!dfind] == "")) stop(
                "When DATA is a list, all its elements that are not data frames must be named."
                )
        # Checking for duplicate names....
        dupnames <- unique(dlnames[duplicated(dlnames)])
        if(length(dupnames) > 0)
            stop(paste(
                "The following names are used more than once in DATA:",
                paste(dupnames, collapse = ", ")))
        if(any(dfind)) {
            dataold <- DATA
            DATA <- vector("list", 0)
            for(i in seq(along = dataold)) {
                if(dfind[i])
                    DATA <- c(DATA, as.list(dataold[[i]]))
                else DATA <- c(DATA, dataold[i])
            }
            dataold <- NULL
        }
        dlnames <- names(DATA)
        dupnames <- unique(dlnames[duplicated(dlnames)])
        # Checking for duplicated names again (now that columns of data frames are included)....
        if(length(dupnames) > 0) stop(paste(
                "The following names are used more than once in DATA (at least once within a data frame):",
                paste(dupnames, collapse = ", ")))
        # Checking for factors....
        factorind <- sapply(DATA, is.factor)
        if(any(factorind))
            stop(paste(
                "DATA may not include factors. One or more factor variables were detected:",
                paste(dlnames[factorind], collapse = ", ")))
        # Checking for character vectors....
        charind <- sapply(DATA, is.character)
        if(any(charind))
            stop(paste(
                "WinBUGS does not handle character data. One or more character variables were detected:",
                paste(dlnames[charind], collapse = ", ")))
        # Checking for complex vectors....
        complexind <- sapply(DATA, is.complex)
        if(any(complexind))
            stop(paste(
                "WinBUGS does not handle complex data. One or more complex variables were detected:",
                paste(dlnames[complexind], collapse = ", ")))
        # Checking for values farther from zero than 1E+38 (which is limit of single precision)....
        toobigind <- sapply(DATA, function(x)
        {
            y <- abs(x[!is.na(x)])
            any(y[y > 0] > 9.9999999999999998e+37)
        }
        )
        if(any(toobigind))
            stop(paste(
                "WinBUGS works in single precision. The following variables contain data outside the range +/-1.0E+38: ",
                paste(dlnames[toobigind], collapse = ", "),
                ".\n", sep = ""))
        # Checking for values in range +/-1.0E-38 (which is limit of single precision)....
        toosmallind <- sapply(DATA, function(x)
        {
            y <- abs(x[!is.na(x)])
            any(y[y > 0] < 9.9999999999999996e-39)
        }
        )
        n <- length(dlnames)
        data.string <- as.list(rep(NA, n))
        for(i in 1:n) {
            if(length(DATA[[i]]) == 1) {
                ac <- toSingleS4(DATA[[i]])
                data.string[[i]] <- paste(names(DATA)[i], "=",
                    ac, sep = "")
                next
            }
            if(is.vector(DATA[[i]]) & length(DATA[[i]]) > 1) {
                ac <- toSingleS4(DATA[[i]])
                data.string[[i]] <- paste(names(DATA)[i], "=c(",
                    paste(ac, collapse = ", "), ")", sep = 
                    "")
                next
            }
            if(is.array(DATA[[i]])) {
                ac <- toSingleS4(aperm(DATA[[i]]))
                data.string[[i]] <- paste(names(DATA)[i], 
                    "= structure(.Data= c(", paste(ac,
                    collapse = ", "), "), \n   .Dim=c(",
                    paste(as.character(dim(DATA[[i]])),
                    collapse = ", "), "))", sep = "")
            }
        }
        data.tofile <- paste("list(", paste(unlist(data.string), 
            collapse = ", "), ")", sep = "")
        if(any(toosmallind))
            warning(paste(
                "WinBUGS works in single precision. The following variables contained nonzero data",
                "\ninside the range +/-1.0E-38 that were set to zero: ",
                paste(dlnames[toosmallind], collapse = ", "),
                ".\n", sep = ""))
        return(data.tofile)
    }
    rslt <- formatDataS4(DATA)
    cat(rslt, file = towhere, fill = fill)
    invisible(0)
}


toSingleS4 <-
#
# Takes numeric vector and removes digit of exponent in scientific notation (if any)
# 
# Written by Terry Elrod. Disclaimer: This function is used at the user's own risk. 
# Please send comments to Terry.Elrod@UAlberta.ca.
# Revision history: 2002-11-19. Fixed to handle missing values properly.
function(x)
{
    xdim <- dim(x)
    x <- as.character(as.single(x))

    # First to look for positives:
    pplus <- regMatchPos(x, "e\\+0")
    pplusind <- apply(pplus, 1, function(y)
    (!any(is.na(y))))
    if(any(pplusind)) {
        # Making sure that periods are in mantissa...
        init <- substring(x[pplusind], 1, pplus[
            pplusind, 1] - 1)
        #...preceeding exponent
        pper <- regMatchPos(init, "\\.")
        pperind <- apply(pper, 1, function(y)
        (all(is.na(y))))
        if(any(pperind))
            init[pperind] <- paste(init[pperind],
                ".0", sep = "")
        # Changing the format of the exponent...
        x[pplusind] <- paste(init, "E+", substring(
            x[pplusind], pplus[pplusind, 2] + 1),
            sep = "")
    }
    # Then to look for negatives:
    pminus <- regMatchPos(x, "e\\-0")
    pminusind <- apply(pminus, 1, function(y)
    (!any(is.na(y))))
    if(any(pminusind)) {
        # Making sure that periods are in mantissa...
        init <- substring(x[pminusind], 1, pminus[
            pminusind, 1] - 1)
        #...preceeding exponent
        pper <- regMatchPos(init, "\\.")
        pperind <- apply(pper, 1, function(y)
        (all(is.na(y))))
        if(any(pperind))
            init[pperind] <- paste(init[pperind],
                ".0", sep = "")
        # Changing the format of the exponent...
        x[pminusind] <- paste(init, "E-", substring(
            x[pminusind], pminus[pminusind, 2] +
            1), sep = "")
    }
    x
}

}
