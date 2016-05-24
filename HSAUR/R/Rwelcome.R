
Rwelcome <- function() {

    tversion <- paste(version$major, version$minor, sep = ".")
    tdate <- paste(version$year, version$month, version$day, sep = "-") 
    x <- c(paste("R : Copyright", version$year, "The R Foundation for Statistical Computing"),
           paste("Version", tversion, paste("(", tdate, "),", sep = ""), 
                 "ISBN 3-900051-07-0"),
           " ",
           "R is free software and comes with ABSOLUTELY NO WARRANTY.",
           "You are welcome to redistribute it under certain conditions.",
           "Type 'license()' or 'licence()' for distribution details.",
           " ",
           "R is a collaborative project with many contributors.",
           "Type 'contributors()' for more information and",
           "'citation()' on how to cite R or R packages in publications.",
           " ",
           "Type 'demo()' for some demos, 'help()' for on-line help, or",
           "'help.start()' for an HTML browser interface to help.",
           "Type 'q()' to quit R.\n")
    cat(paste(x, collapse = "\n"))
}

exename <- function() {

    tversion <- paste(version$major, "0", substr(version$minor, 1, 1),
                      substr(version$minor,3,3), sep = "")
    return(paste("rw", tversion, ".exe", sep = ""))
}
