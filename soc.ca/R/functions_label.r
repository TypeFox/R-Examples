## Labels 

#' Add values to label
#' 
#' Adds values to the end of the label of each modality.
#' @param  object is a soc.ca object
#' @param  value the type of values added to the labels. "freq" adds
#'   frequencies, "mass" adds mass values to the active modalities, "ctr" adds contribution values to the active modalities, "cor" adds correlation values.
#'   value also accepts any vector with the length of the number of active
#'   modalities. "linebreak" adds a linebreak \code{\\n} after the first ":" in the label.
#' @param  prefix if "default" an appropriate prefix is used
#' @param  suffix the suffix
#' @param  dim the dimension from which values are retrieved
#' @return a soc.ca object with altered labels in names.mod and names.sup
#' @export
#' @examples
#' example(soc.ca)
#' result.label  <- add.to.label(result)
#' result.label$names.mod
#' result.label  <- add.to.label(result, value = "ctr", dim = 2)
#' result.label$names.mod
#' result.label  <- add.to.label(result, value = result$variable, prefix = " - ", suffix = "")
#' result.label$names.mod
#' result.label  <- add.to.label(result, value = "linebreak")
#' result.label$names.mod
#' map.ctr(result.label)

add.to.label <- function(object, value = "freq", prefix = "default", suffix = ")", dim = 1){
  
  # Names
  names.mod           <- object$names.mod
  names.sup           <- object$names.sup
  
  # Prefix
  if (identical(prefix, "default") & identical(value, "freq")) prefix   <- " (n:"
  if (identical(prefix, "default") & identical(value, "mass")) prefix   <- " (mass:"
  if (identical(prefix, "default") & identical(value, "ctr")) prefix    <- " (ctr:"
  if (identical(prefix, "default") & identical(value, "cor")) prefix    <- " (cor:"
  if (identical(prefix, "default") & length(value) > 1 )       prefix   <- " ("
  
  # Values
  if (identical(value, "freq")){
  val.mod             <- object$freq.mod
  val.sup             <- object$freq.sup
  object$names.mod    <- paste(names.mod, prefix, val.mod, suffix, sep = "")
  object$names.sup    <- paste(names.sup, prefix, val.sup, suffix, sep = "")
  }
  
  if (identical(value, "ctr")){
  val.mod             <- round(object$ctr.mod[dim, ] * 100, 1)
  object$names.mod    <- paste(names.mod, prefix, val.mod, suffix, sep = "")
  }
  
  if (identical(value, "mass")){
    val.mod           <- round(object$mass.mod * 100, 1)
    object$names.mod  <- paste(names.mod, prefix, val.mod, suffix, sep = "")
  }
  
  if (identical(value, "cor")){
    val.mod           <- round(object$cor.mod[, dim], 2)
    val.sup           <- round(object$cor.sup[, dim], 2)
    object$names.mod  <- paste(names.mod, prefix, val.mod, suffix, sep = "")
    object$names.sup  <- paste(names.sup, prefix, val.sup, suffix, sep = "")
  }
  
  if (length(value) > 1){
    object$names.mod  <- paste(names.mod, prefix, value, suffix, sep = "")
  }

  if (identical(value, "linebreak")){
    object$names.mod  <- sub(":", ":\n", names.mod)
    object$names.sup  <- sub(":", ":\n", names.sup)
    
  }
  
  colnames(object$indicator.matrix) <- object$names.mod
  
  return(object)
}

#' Exports the labels of a soc.ca object into a csv file.
#' 
#' This function allows easy translation and renaming of modalities by exporting
#' the labels into a .csv file that is easier to work with.
#' 
#' Two columns are created within the .csv: 'New label' and 'Old label'. In the
#' 'New label' column you write the new labels. Remember to leave 'Old label'
#' unchanged as this column is used for matching.
#' 
#' If you want to add frequencies to the labels with the \link{add.to.label}
#' function you should do this after exporting and assigning labels with the
#' \link{assign.label} function. Otherwise the matching of the labels is likely
#' to fail.
#' @param object is a soc.ca object
#' @param file is the name and path of the exported file
#' @param encoding is the character encoding of the exported file
#' @param overwrite decides whether to overwrite already existing files
#' @return A .csv with two columns and preferably UTF-8 encoding.
#' @export

export.label    <- function(object, file = FALSE, encoding = "UTF-8", overwrite = FALSE){
  
  names         <- c(object$names.mod, object$names.sup, object$names.ind)
  ca.label      <- cbind(names, names)
  colnames(ca.label)  <- c("New.label", "Old.label")
  
  if (identical(file, FALSE) == TRUE){
    file    <- paste("label_",deparse(substitute(object)), ".csv", sep = "")
  }
  if(file.exists(file) == TRUE & identical(overwrite, FALSE)) stop("File already exists")
  write.csv(ca.label, file = file, fileEncoding = encoding)

}


#' Assign new labels
#' 
#' Assigns new labels to a soc.ca object. The input labels are defined in a .csv
#' file created by the \link{export.label} function.
#' @param object is a soc.ca object
#' @param file is the path of the .csv file with the new labels. The file is
#'   preferably created by the \link{export.label} function
#' @param encoding is the encoding of the imported file
#' @param sep is the seperator used to create the imported .csv file
#' @return a soc.ca object with altered labels in \code{object$names.mod}, \code{object$names.ind} and
#'   \code{object$names.sup}
#' @details To use this function first export the labels from your soc.mca
#'   analysis with the \link{export.label} function. Then open and edit the
#'   created file with your favorite spreadsheet editor, like LibreOffice Calc.
#'   Change labels in the "new.label" column to the desired values and save. Use the
#'   assign.label function but remember to assign the results into a new object
#'   or overwrite the existing object.
#' @seealso \link{export.label}, \link{add.to.label}
#' @export

assign.label <- function(object, file = FALSE, encoding = "UTF-8", sep = ","){
  if (identical(file, FALSE) == TRUE){
    file <- paste("label_", deparse(substitute(object)), ".csv", sep = "")
  }
  label     <- read.csv(file, encoding = encoding, sep  =  sep)
  
  names.mod <- as.character(object$names.mod)
  names.sup <- as.character(object$names.sup)
  names.ind <- as.character(object$names.ind)
  
  new.label <- as.character(label$New.label)
  old.label <- as.character(label$Old.label)
  
  for (i in 1:length(old.label)){
    indices   <- which(old.label[i] == names.mod)
    names.mod[indices] <- new.label[i]
  }
  
  for (i in 1:length(old.label)){
    indices   <- which(old.label[i] == names.sup)
    names.sup[indices] <- new.label[i]
  }
  
  for (i in 1:length(old.label)){
    indices   <- which(old.label[i] == names.ind)
    names.ind[indices] <- new.label[i]
  }
  
  
  object$names.mod <- names.mod
  object$names.sup <- names.sup
  object$names.ind <- names.ind
  return(object)

}