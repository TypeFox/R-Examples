fences.summary <-
function (group, x, file = NULL, units = "ppm") 
{
    # Function to generate files of framework summary statistics to be saved,
    # by default in the R working directory.  If provided "file" should name 
    # the folder where the results are to be saved.  The function appends the 
    # data frame, framework grouping and variable names together with ".csv" 
    # prior to opening a file for the results; these can be opened/viewed 
    # with MS Excel.  To function correctly the data frame containing the
    # data must be attached prior to running this function, if it has not
    # already been done so, i.e. attach(dfname), and at close it should be
    # detached if necesary, detach(dfname).
    #
    # To function correctly the data frame containing the data must be attached 
    # prior to running this function, if it has not already been done so, i.e. 
    # attach(dfname), and at close it should be detached if necesary, detach(dfname).
    #
    # NOTE: Prior to using this function the data frame/matrix containing the
    # variable, 'x', data must be run through ltdl.fix.df to convert any <dl
    # -ve values to positive half that value, and set zero2na = TRUE if it is
    # required, to convert any zero values or other numeric codes representing 
    # blanks to NAs.
    #
    # Functions fences and summary.stats is used to compute the fences. 
    #
    dfname <- search()[[2]]
    groupname <- deparse(substitute(group))
    xname <- deparse(substitute(x))
    filename <- paste(dfname, "_", groupname, "_", xname, "_fences.txt", sep = "")

    if (is.null(file)) {
        wdname <- getwd()
        cat("  By default the output file will be saved in the R working directory as:\n",
            "  \"dataframename_groupname_xname_fences.txt\"\n",
            " To specify an alternate location and file name prefix, to which:\n",
            "  \"/dfname_groupname_xname_fences.txt\" will be appended, set, for example:\n",
            "  file = \"D://R_work//Project3//C_soils\"\n")
        filename <- paste(wdname, "/", filename, sep = "")
    }
    else filename <- paste(file, "/", filename, sep = "")  
    #
    cat("\n  Fences for", xname, "subset by", groupname, "- output will be in:\n  ", 
        filename, "\n\n")
    sink(filename)
    on.exit(sink())
    #
    framework.fences <- tapply(x, group, fences, units = units, 
        display = FALSE)
    nfences <- length(framework.fences)
    for (i in 1:nfences) {
        gi <- names(framework.fences[i])
        titl <- paste(xname, "[", gi, "]", " (Units = ", units, ")", sep = "")
        table <- unlist(framework.fences[i], use.names = FALSE)
        cat(titl, "\tN =", table[1], "\tNAs =", table[2], 
            "\t2%ile =", table[14], "\t98%ile =", table[10], 
            "\n\t Mean\t SD\t Median\t MAD\t\t Mean\t Med\t Tukey Fences\n",
            "\t\t\t\t\t\t \2612SD\t \2612MAD\t\t(actual)\n", 
            "\t", table[3], "\t", table[4], "\t", table[5], "\t", table[6], "\t\t",
            table[7], "\t", table[8], "\t", table[9], "\t", table[25], "\n\t\t\t\t\t\t", 
            table[11], "\t", table[12], "\t", table[13], "\t", table[26], 
            "\nLog10\t", table[15], "\t", table[16], "\t", table[17], "\t", table[18], "\t\t", 
            table[19], "\t", table[20], "\t", table[21], "\t", table[27], "\n\t\t\t\t\t\t",
            table[22], "\t", table[23], "\t", table[24], "\t", table[28], 
            "\nLogit\t", table[29], "\t", table[30], "\t", table[31], "\t", table[32], "\t\t",
            table[33], "\t", table[34], "\t", table[35], "\t", table[39], "\n\t\t\t\t\t\t", 
            table[36], "\t", table[37], "\t", table[38], "\t", table[40], "\n\n")
    }
    invisible()
}
