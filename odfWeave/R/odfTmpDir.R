# This generates the path of a directory to use
# in odfWeave.  It won't exist, as it is the job
# of the odfWeave function to create it.
odfTmpDir <- function()
{
   tmpPath <- tempdir()
   suffix <- paste(
      "odfWeave", 
      format(Sys.time(), "%d%H%M%S"), 
      round(runif(1)*1000, 0),
      sep = "")
   paste(tmpPath, "/", suffix, sep = "")
}

