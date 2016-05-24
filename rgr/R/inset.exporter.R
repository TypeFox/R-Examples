inset.exporter <-
function(x, xlab = deparse(substitute(x)), log = FALSE, xlim = NULL, nclass = NULL, 
     ifnright = TRUE, file = NULL, table.cex = 0.7, gtype = "emf" , ...)
{
     # Wrapper to save the graphics output from a run of inset.
     #
     # NOTE: This functionality is only available under the Windows OS.
     #
     # To function correctly the data frame containing the data must be attached 
     # prior to running this function, if it has not already been done so, i.e. 
     # attach(dfname), and at close it should be detached if necesary, detach(dfname).
     #
     # NOTE: Prior to using this function the data frame/matrix containing the
     # variable, 'x', data must be run through ltdl.fix.df to convert any <dl
     # -ve values to positive half that value, and set zero2na = TRUE if it is
     # required to convert any zero values or other numeric codes representing 
     # blanks to NAs.
     #
     info <- Sys.info()
     if(!info[1] == "Windows") stop("\n  This function only available under the Windows OS\n")
     ##
     dfname <- search()[[2]]
     xname <- deparse(substitute(x))
     if (is.null(file)) {
         file <- getwd()
         cat("  By default the output file will be saved in the R working directory as:\n",
             "  \"dataframename_inset_xname.emf\"\n",
             " To specify an alternate location and file name prefix, to which:\n",
             "  \"dfname_inset_xname.emf\" will be appended, set, for example:\n",
             "  file = \"D://R_work//Project3\"\n\n")
     }
     file.name <- paste(file, "/", dfname, "_inset_", xname, ".", gtype, sep = "")
     cat("  Graphics file will be saved as:", file.name, "\n\n")
     ##
     inset(x, xlab = xlab, log = log, xlim = xlim, nclass = nclass, ifnright = ifnright,
         table.cex = table.cex, ...)
     savePlot(filename = file.name, type = gtype)
     invisible()
}
