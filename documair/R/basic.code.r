
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
analyse8description <- function(dfile)
#TITLE returns useful component from a description file 
#DESCRIPTION \samp{dfile} must be a valid \samp{DESCRIPTION}
# file of a package. It is read and some useful information
# about the package are returned.
#DETAILS
# Not sure that the analysis be completely general. 
#KEYWORDS description extract
#INPUTS
#{dfile}<< A \samp{character(1)} indicating the file to explore.>>
#[INPUTS]
#VALUE
# A named list of version, date and dependencies.
#EXAMPLE
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE
#AUTHOR J.-B. Denis
#CREATED 14_05_22
#REVISED 14_05_22
#--------------------------------------------
{
  # constant
  champs <- c("Package","Version","Date",
              "Title","Author","Depends");
  # checking
  rrrobject9(dfile,"character",1,
             mensaje="'dfile' must be a 'character(1)'");
  rrrobject9(dfile,"file",1,
             mensaje="'dfile' must designate a file");
  # reading it
  descri <- rrrfile2text(dfile);
  # initialization
  res <- vector("list",0);
  # getting something from it
  for (lig in descri) {
    decode <- strsplit(lig,":",fixed=TRUE)[[1]];
    clef <- decode[1];
    if (nchar(clef) > 0) {
      serr <- which(clef==champs);
      if (length(serr)==1) {
        res[[clef]] <- rrrform3crop(decode[2]);
      }
    }
  }
  # returning
  res;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
