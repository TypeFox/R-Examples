###
### $Id: fullfile.R 48 2014-02-05 20:50:54Z plroebuck $
###
### Build full filename from parts.
###


##-----------------------------------------------------------------------------
fullfile <- function(...) {
    return(file.path(...))
}

