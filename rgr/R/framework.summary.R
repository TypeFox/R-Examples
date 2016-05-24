framework.summary <-
function(group, x, file = NULL)
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
     # NOTE: Prior to using this function the data frame/matrix containing the
     # variable, 'x', data must be run through ltdl.fix.df to convert any <dl
     # -ve values to positive half that value, and set zero2na = TRUE if it is
     # required, to convert any zero values or other numeric codes representing 
     # blanks to NAs.
     #
     # Function framework.stats and summary.stats are used to compute the 
     # summary statistics. 
     #
     dfname <- search()[[2]]
     groupname <- deparse(substitute(group))
     xname <- deparse(substitute(x))
     filename <- paste(dfname, "_", groupname, "_", xname, ".csv", sep = "")
     if (is.null(file)) {
         file <- getwd()
         cat("  By default the output file will be saved in the R working directory as:\n",
             "  \"dataframename_groupname_xname.csv\"\n",
             " To specify an alternate location and file name prefix, to which:\n",
             "  \"/dataframename_groupname_xname.csv\" will be appended, set, for example:\n",
             "  file = \"D://R_work//Project3\"\n")
     }
     filename <- paste(file, "/", filename, sep = "")  
     #
     cat("\n  Summary statistics for variable", xname, "subset by", groupname, "- output will be in:\n  ", 
         filename, "\n\n")
     sink(filename)
     on.exit(sink())
     #
     cat("Variable,Group,N,NA,Min,2%ile,5%ile,10%ile,25%ile,Median,75%ile,90%ile,95%ile,98%ile,Max,LCI,UCI,MAD,IQSD,Mean,SD,CV%"
         )
     #
     framework.stats <- tapply(x, group, framework.stats)
     nstats <- length(framework.stats)
     for(i in 1:nstats) {
         gi <- names(framework.stats[i])
         ii <- unlist(framework.stats[i], use.names = FALSE)
         cat("\n", xname, ",", gi, ",", ii[1], ",", ii[2], ",", ii[3], ",", ii[4], ",", ii[5], 
          ",", ii[6],",", ii[7], ",", ii[8], ",", ii[9], ",", ii[10], ",", ii[11], ",", ii[12],
          ",", ii[13], ",", ii[14], ",", ii[15], ",", ii[16], ",", ii[17], ",", ii[18], ",",
          ii[19], ",", ii[20], sep = "")
     }
     invisible()
}
