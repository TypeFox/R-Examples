###############################################
# --------------------------------------------#
# saveRand                                    #
# --------------------------------------------#
###############################################

# --------------------------------------------
# Function for saving the parameters of an randSeq object
# --------------------------------------------

#' Saving a randomization lists
#' 
#' Saves the parameters of a \code{randSeq} object in a \code{.csv} data
#' sheet.
#' 
#' @family saving functions
#' 
#' @inheritParams overview
#' @param obj object of class \code{randSeq} containing a single randomization sequence.
#'
#' @return Creates a \code{.csv} data in the home folder and saves the randomization list as a column vector.
#'
#' @export
saveRand <- function(obj, file = "randList.csv") {
  if(!("randSeq" %in% is(obj))) stop("Object not of class randSeq")
  if (nrow(obj@M) > 1) stop("Object should only contain one randomization sequence.")
  
  write(paste("This document was generated on",  format(Sys.time(), "%a %b %d %Y"),
              "at",  format(Sys.time(), "%X") ,".\n"),
        file = file)

  sessInfo <- sessionInfo()
  write(paste("The randomizeR package of version",  packageVersion("randomizeR"),
              "was used for generating the randomization list with the",
              sessInfo$R.version$version.string,".") ,
          file = file, append = TRUE)
  
  write(paste("\nRandomization Method:", getDesign(obj), "\n"),
        file = file, append = TRUE)
  # iterate through all slots of the object
  names <- slotNames(obj)
  # without M
  names <- names[!(names == "M")] 
  for(name in names) {
    write(paste(name, ":\t", paste(slot(obj, name), collapse = ", "), sep=""),
          file = file, append = TRUE)
  }
  
  write(paste("\nLegend: \nN := number of included patients ",
        "\nK := number of treatment groups",
        "\nratio := allocation ratio of the trial",
        "\ngroups := names for the investigated groups"),
         file = file, append = TRUE)
  write("For specific randomization parameters see the help of the randomizeR package.",
         file = file, append = TRUE)

  write("\n" , file = file, append = TRUE)
  randList <- getRandList(obj)
  
  dataFrame <- data.frame("Allocation" = 1:length(randList))
  dataFrame$Treatment = as.vector(randList)
  
  suppressWarnings(write.table(dataFrame, file = file,  sep = "\t", append = TRUE, 
                               row.names = FALSE, col.names = TRUE))
  
  # old version (output is transposed)
  #colnames(randList) <- paste("Allocation", 1:ncol(randList))
  #write(colnames(randList) , file = file,  sep = "\t", ncolumns = ncol(randList),
  #      append = TRUE)
  #write(t(randList) , file = file,  sep = "\t", ncolumns = ncol(randList),
  #      append = TRUE)
  
}











