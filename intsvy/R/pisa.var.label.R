pisa.var.label <- 
function(folder=getwd(), student.file, parent.file=c(), school.file=c(), name="Variable labels", output=getwd()) {
  
  # Student file required for country labels
  if(missing(student.file)) {
    stop("the student file is required")
  }
  
  if(!missing(folder)) {
    student.file =  file.path(folder, paste(student.file, sep=""))
  }
    
  if(!missing(folder) & !missing(school.file)) {
    school.file =   file.path(folder, paste(school.file, sep=""))
  }

  if(!missing(folder) & !missing(parent.file)) {
    parent.file =  file.path(folder, paste(parent.file, sep=""))
  }
  
  # Retrieve file name
  files.all <- list(student.file, parent.file, school.file)
  names(files.all) <- c('Student', 'Parent', 'School')
  
  # Remove null elements in list
  files.all <- files.all[lapply(files.all, length)>0]
  
  # Retrieve var labels
  var.label <- lapply(files.all, function(x) description(spss.system.file(x[[1]], to.lower=FALSE)))  
  
  # Read student file and participating countries 
  country <- names(table(spss.system.file(files.all[["Student"]], to.lower=FALSE)[,"CNT"]))
  
  # Participating countries in dataset (must be all)
  country.list <- pisa.country[pisa.country[, "ISO"] %in% country, ]
  rownames(country.list) <-NULL
  
  # setdiff(country[,1], pisa.country$ISO) must be zero
  
  var.label[[length(files.all)+1]] <- country.list
  names(var.label)[length(var.label)] <-"Participating countries"
  
  # Print labels in list and text file
  capture.output(var.label, file=file.path(output, paste(name, ".txt", sep="")))
  cat('The file "', paste(name, ".txt", sep=""), '" in directory "', output, '" contains the variable labels of the complete dataset', sep=' ', "\n")
  return(var.label)
}
