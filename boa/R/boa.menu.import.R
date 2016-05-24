"boa.menu.import" <-
function()
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
{
   mtitle <- "\nIMPORT DATA MENU\n----------------"
   choices <- c("Back",
                "---------------------------+",
                "CODA Output Files          |",
                "Flat ASCII Files           |",
                "Data Matrix Objects        |",
                "View Format Specifications |",
                "Options...                 |",
                "---------------------------+")
   idx <- 1
   while(idx > 0) {
      idx <- menu(choices, title = mtitle)
      switch(idx,
         "1" = idx <- -1,
         "2" = NULL,
         "3" = { cat("\nEnter index filename prefix without the .ind ",
                     "extension [Working Directory: ",
                     deparse(boa.par("path")), "]\n", sep = "")
                 value <- scan(what = "", n = 1, strip.white = TRUE)
                 cat("\nEnter output filename prefix without the .out ",
                     "extension [Default: \"", value, "\"]\n", sep = "")
                 value <- c(value, scan(what = "", n = 1, strip.white = TRUE))
                 if(boa.chain.import(value, type = "BUGS"))
                    cat("+++ Data successfully imported +++\n")
               },
         "4" = { cat("\nEnter filename prefix without the ",
                     boa.par("ASCIIext"), " extension [Working Directory: ",
                     deparse(boa.par("path")), "]\n", sep = "")
                 value <- scan(what = "", n = 1, strip.white = TRUE)
                 if(boa.chain.import(value, type = "ASCII"))
                    cat("+++ Data successfully imported +++\n")
               },
         "5" = { cat("\nEnter object name [none]\n")
                 value <- scan(what = "", n = 1, strip.white = TRUE)
                 if(boa.chain.import(value, type = "S"))
                    cat("+++ Object successfully imported +++\n")
               },
         "6" = { cat("\nCODA\n",
                     "- CODA index (*.ind) and output (*.out) files produced",
                     " by WinBUGS\n",
                     "- Index and output files must be saved as ASCII text\n",
                     "- Files must be located in the Working Directory (see",
                     " Options)\n",
                     "\nASCII\n",
                     "- ASCII file (*", boa.par("ASCIIext"), ") containing",
                     " the monitored parameters from one run of the\n",
                     "  sampler\n",
                     "- Parameters are stored in space or tab delimited",
                     " columns\n",
                     "- Parameter names must appear in the first row\n",
                     "- Iteration numbers may be specified in a column",
                     " labeled 'iter'\n",
                     "- File must be located in the Working Directory (see",
                     " Options)\n",
                     "\nMatrix Object\n",
                     "- R numeric matrix whose columns contain the monitored",
                     " parameters from one\n",
                     "  run of the sampler\n",
                     "- Iteration numbers and parameter names may be",
                     " specified in the dimnames\n", sep = "")
                 cat("\nPress <ENTER> to continue")
                 readline()
               },
         "7" = boa.menu.setpar("Import"),
         "8" = NULL
      )
   }

   return(abs(idx))
}

