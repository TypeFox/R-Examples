###############################################
# --------------------------------------------#
# saveAssess                                 #
# --------------------------------------------#
###############################################

# --------------------------------------------
# Function for saving the parameters of an Asssess object
# --------------------------------------------

#' Saving an assess object
#' 
#' Saves the full information of an \code{assess} object in a \code{.csv} data
#' sheet.
#' 
#' @family saving functions
#' 
#' @inheritParams overview
#' @param obj object of class assessmentment, e.g. the output of the \code{\link{assess}} function.
#'
#' @return Creates a \code{.csv} data in the home folder.
#'
#' @export
saveAssess <- function(obj, file = "assessObject.csv") {
  if(!("assessment" %in% is(obj))) stop("Object not of class asssessment")
  
  write(paste("This document was generated on",  format(Sys.time(), "%a %b %d %Y"),
              "at",  format(Sys.time(), "%X") ,".\n"),
        file = file)
  
  sessInfo <- sessionInfo()
  write(paste("The randomizeR package of version",  packageVersion("randomizeR"),
              "was used for generating this object with the",
              sessInfo$R.version$version.string,".") ,
        file = file, append = TRUE)
  
  write(paste("\nRandomization Method:", obj@design, "\n"),
        file = file, append = TRUE)
  # iterate through all slots of the object
  names <- slotNames(obj)
  names <- names[!(names %in% c("D", "design"))] 
  for(name in names) {
    write(paste(name, ":\t", paste(slot(obj, name), collapse = ", "), sep=""),
          file = file, append = TRUE)
  }
  
  write(paste("\nLegend: \nN := number of included patients ",
              "\nK := number of treatment groups",
              "\ngroups := names for the investigated groups"),
        file = file, append = TRUE)
  write("For specific randomization parameters see the help of the randomizeR package.",
        file = file, append = TRUE)
  
  write("\n" , file = file, append = TRUE)
  write(colnames(obj@D) , file = file,  sep = "\t", ncolumns = ncol(obj@D),
        append = TRUE)
  
  obj@D[,-1] <- round(obj@D[,-1], digits = 3)

  write.table(obj@D, file = file, sep = "\t", append = TRUE, row.names = FALSE, col.names = FALSE)
  
}
