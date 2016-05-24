### Tools ######################################################################

##' @title Check whether a System Command Exists
##' @param command character string containing the command
##' @return logical
##' @author Marius Hofert
##' @note See http://stackoverflow.com/questions/592620/check-if-a-program-exists-from-a-bash-script
command.exists <- function(command)
    system(paste("command -v", command, ">/dev/null 2>&1")) == 0

##' @title Returning the Extension of a File
##' @param x A file name or path
##' @return The file extension
##' @author Marius Hofert
##' @note Adapted from tools::file_ext()
file_ext <- function(x)
{
    pos <- regexpr("\\.([[:alnum:]]+)$", x)
    ifelse(pos > -1L, substring(x, pos + 1L), "")
}

##' @title Returning the Basenmane without Extension of a File
##' @param x A file name or path
##' @return The file basename without extension
##' @author Marius Hofert
##' @note Adapted from tools::file_path_sans_ext()
file_path_sans_ext <- function (x)
    sub("([^.]+)\\.[[:alnum:]]+$", "\\1", x)

##' @title Add a String to the End of a File Name (before the Extension)
##' @param x A file name or path
##' @param s A string to be added
##' @return New file name including s
##' @author Marius Hofert
add_to_name <- function(x, s)
    paste0(file_path_sans_ext(x),s,".",file_ext(x))

