R.to.Tcl <-
function (character.vector)
##  converts a character vector into a brace-delimited Tcl list
{
    if (length(character.vector) == 0) list()
    else paste("{", paste(character.vector, collapse = "} {"),
               "}",
               sep = "")
}

Tcl.to.R <-
function (tcl.list)
##  converts a fully brace-delimited Tcl list into a character
##  vector in R
{
    tcl.list <- substring(tcl.list, 2, nchar(tcl.list) - 1)
    strsplit(tcl.list, split = "} {", fixed = TRUE)[[1]]
}

