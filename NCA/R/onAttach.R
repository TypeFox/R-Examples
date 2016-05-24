.onAttach <- function(...) {

    meta  <- packageDescription("NCA")
    year  <- sub("-.*", "", meta$Date)
    title <- meta$Title

    msg1  <- sprintf("Dul, J. %s.", year)
    msg2  <- sprintf("%s.", title)
    msg3  <- sprintf("R Package Version %s.\n", meta$Version)
    msg4  <- "URL: http://cran.r-project.org/web/packages/NCA/"

    msg5  <- "This package is based on:"
    msg6  <- "Dul, J. (2015) Logic and Methodology of"
    msg7  <- "'Necessary but not Sufficient' causality."
    msg8  <- "Published in Organizational Research Methods (Sage publishers)"
    msg9  <- "http://dx.doi.org/10.2139/ssrn.2588480"

    msg10  <- "A BibTeX entry is provided by:"
    msg11  <- 'citation("NCA")'

    packageStartupMessage("\nPlease cite the NCA package as:\n\n",
                          strwrap(msg1, indent = 2, exdent = 2), "\n",
                          strwrap(msg2, indent = 2, exdent = 2), "\n",
                          strwrap(msg3, indent = 2, exdent = 2), "\n",
                          strwrap(msg4, indent = 2, exdent = 2), "\n",
                          "\n",
                          strwrap(msg5, indent = 2, exdent = 2), "\n",
                          strwrap(msg6, indent = 2, exdent = 2), "\n",
                          strwrap(msg7, indent = 2, exdent = 2), "\n",
                          strwrap(msg8, indent = 2, exdent = 2), "\n",
                          strwrap(msg9, indent = 2, exdent = 2), "\n",
                          "\n",
                          strwrap(msg10, indent = 0, exdent = 2), "\n",
                          strwrap(msg11, indent = 2, exdent = 2), "\n",
                          "\n")
}

