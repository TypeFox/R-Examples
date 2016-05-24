pisa.select.merge <- 
function(folder=getwd(), student.file, parent.file=c(), school.file=c(), countries, student=c(), parent, school) {
  
  # Remove leading and trailing whitespaces in var labes  
  if(!missing(student) & !is.null(student)) {
    student = gsub("^[[:space:]]+|[[:space:]]+$", "", student)
    student =toupper(student)
  }
  
  if(!missing(parent)){
    parent = gsub("^[[:space:]]+|[[:space:]]+$", "", parent)
    parent = toupper(parent)
  }
  
  if(!missing(school)){
    school = gsub("^[[:space:]]+|[[:space:]]+$", "", school)
    school = toupper(school)
  }
  
  # No variables selected (error)  
  
  if (missing(student) & missing(parent) & missing(school)) {
    stop("no variables are selected")
  }
  
  # Creating student, parent, and school filepaths 
  
  # Student file required for country labels
  if(missing(student.file)) {
  stop("the student file is required")
  }
  
  if(!missing(folder) & !missing(student.file)) {
    student.file =   file.path(folder, paste(student.file, sep=""))
  }
  
  if(!missing(folder) & !missing(school.file)) {
    school.file =  file.path(folder, paste(school.file, sep=""))
  }
  
  if(!missing(folder) & !missing(parent.file)) {
    parent.file =  file.path(folder, paste(parent.file, sep=""))
  }
  
  # Retrieve file names
  files.all <- list(student.file, parent.file, school.file)
  names(files.all) <- c('Student', 'Parent', 'School')
  
  # Remove null elements in list
  files.all <- files.all[lapply(files.all, length)>0]
  
  # Participating countries (from student file)
  pisa.student <- spss.system.file(files.all[["Student"]], to.lower=FALSE)
  country <- names(table(pisa.student[,"CNT"]))
  
  # If countries missing, all countries selected
  if (missing(countries)) {
    countries <- country
  }
  
  # If countries are entered numerically, change to ISO labels for file selection (next)
  if (is.numeric(countries)) {
    countries=pisa.country[pisa.country$Code %in% countries, "ISO"]
  }
  
  # Have to use spss.system.file, otherwise read.spss crashes
  
  # Student data (need for school and parent data too)
  
  if (!missing(student) | !missing(parent)) {
    
    names(pisa.student) <- toupper(names(pisa.student)) # because stidstd is lowercase sometimes
    
    student.data <- pisa.student[pisa.student["CNT"] %in% countries, 
                                 c("CNT", "SCHOOLID", "STIDSTD",  
                                   names(pisa.student)[grep("^PV", names(pisa.student))], student,
                                   names(pisa.student)[grep("^W_F", names(pisa.student))])]
    
    student.data <- as.data.frame(adj.measlev(student.data))
    
  }
  
  
  # Parental questionnaire
  
  if (!missing(parent)) {
    
    if (is.null(files.all[["Parent"]])) {
      stop("cannot locate parental questionnaire data file")
    }
    
    
    pisa.parent <- spss.system.file(files.all[["Parent"]], to.lower=FALSE)
    names(pisa.parent) <- toupper(names(pisa.parent))
    
    parent.data <- pisa.parent[pisa.parent["CNT"] %in% countries, 
                               c("CNT", "SCHOOLID", "STIDSTD", parent)]
    
    parent.data <- as.data.frame(adj.measlev(parent.data))
    
  }
  
  
  # School data
  
  if (!missing(school)) {
    
    if (is.null(files.all[["School"]])) {
      stop("cannot locate school data file")
    }
    
    
    pisa.school <- spss.system.file(files.all[["School"]], to.lower=FALSE)
    names(pisa.school) <- toupper(names(pisa.school))
    
    school.data <- pisa.school[pisa.school["CNT"] %in% countries, 
                               c("CNT", "OECD", "SCHOOLID", school, "W_FSCHWT")]
    
    school.data <- as.data.frame(adj.measlev(school.data))
    
  }
  
  
  # Merging data depending on existing datasets/arguments
  
  # Student data available
  
  if (!missing(student) & missing(parent) & missing(school)) {
    pisa.all <- student.data
  }
  
  if (!missing(student) & !missing(parent) & missing(school)) {
    pisa.all <- merge(student.data, parent.data, all.x=TRUE, by=c("CNT", "SCHOOLID", "STIDSTD"))
    pisa.all[, grep("*.y", names(pisa.all))] <- NULL
    names(pisa.all) <- gsub("*.x", "", names(pisa.all))
  }
  
  if (!missing(student) & missing(parent) & !missing(school)) {
    pisa.all <- merge(student.data, school.data, all.x=TRUE, by=c("CNT", "SCHOOLID"))
    pisa.all[, grep("*.y", names(pisa.all))] <- NULL
    names(pisa.all) <- gsub("*.x", "", names(pisa.all))
  }
  
  if (!missing(student) & !missing(parent) & !missing(school)) {
    pisa.all <- merge(student.data, parent.data, all.x=TRUE, by=c("CNT", "SCHOOLID", "STIDSTD"))
    pisa.all <- merge(pisa.all, school.data, all.x=TRUE, by=c("CNT", "SCHOOLID"))
    pisa.all[, grep("*.y", names(pisa.all))] <- NULL
    names(pisa.all) <- gsub("*.x", "", names(pisa.all))
  }
  
  # Parent data available
  
  if (is.null(student) & !missing(parent) & missing(school)) {
    pisa.all <- merge(student.data, parent.data, by=c("CNT", "SCHOOLID", "STIDSTD"))
    pisa.all[, grep("*.y", names(pisa.all))] <- NULL
    names(pisa.all) <- gsub("*.x", "", names(pisa.all))
  }
  
  if (is.null(student) & !missing(parent) & !missing(school)) {
    pisa.all <- merge(school.data, parent.data, by=c("CNT", "SCHOOLID"))
    pisa.all[, grep("*.y", names(pisa.all))] <- NULL
    names(pisa.all) <- gsub("*.x", "", names(pisa.all))
  }
  
  # School data available
  
  if (is.null(student) & missing(parent) & !missing(school)) {
    pisa.all <- school.data
  }
  
  # Create country label variable (not possible to add labels to numeric factor, see to do list)
  pisa.all$IDCNTRYL <- pisa.country[match(pisa.all$CNT, pisa.country$ISO), "Country"]
    
  return(pisa.all)
}
