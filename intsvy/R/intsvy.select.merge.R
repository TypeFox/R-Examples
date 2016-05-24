intsvy.select.merge <-
function(folder=getwd(), countries, student=c(), home, school, teacher, use.labels=TRUE,
         config) {
  
  # Remove leading and trailing whitespaces in var labels  
  if(!missing(student) & !is.null(student)) {
    student = gsub("^[[:space:]]+|[[:space:]]+$", "", student)
  }
  if(!missing(home)){
    home = gsub("^[[:space:]]+|[[:space:]]+$", "", home)
  }
  if(!missing(school)){
    school = gsub("^[[:space:]]+|[[:space:]]+$", "", school)
  }
  if(!missing(teacher)){
    teacher = gsub("^[[:space:]]+|[[:space:]]+$", "", teacher)
  }

  # No variables selected (error)
  if (missing(student) & missing(school) & missing(home) & missing(teacher)) {
    stop("no variables are selected")
  }
  
  # Look for datasets (student, home, school) including student-teacher linkage (.st)
  files.all <- lapply(config$input$prefixes, function(x) list.files(folder, 
  full.names= TRUE, pattern=paste("^", x, ".*.sav$", sep=""), recursive=TRUE))

  # Remove empty elements in list
  files.all <- files.all[lapply(files.all, length) >0]
  
  
  # Name list for identification later, rather than using numbers
  names(files.all) <- unique(unlist(lapply(files.all, function(x) 
    substr(x, nchar(x)+config$input$type_part[1], nchar(x)+config$input$type_part[2]))))
    
  # Country abbrevation labels from existing file names (datasets)
  cntlab <- toupper(unique(unlist(lapply(files.all, function(x) 
    substr(x, nchar(x)+config$input$cnt_part[1], nchar(x)+config$input$cnt_part[2]))))) 
  
  # setdiff(cntlab, iea.country$ISO) needs be zero! all elements in data labels are in userguide
  # Countries in the datasets and userguide
  country.list <- iea.country[iea.country[["ISO"]] %in% intersect(iea.country[["ISO"]], cntlab), ]
  
  # If no countries are selected: select all
  if (missing(countries)) { 
    countries=cntlab
  }
  
  # If countries are entered numerically, change to ISO labels for file selection (next)
  if (is.numeric(countries)) {
    countries=country.list[country.list$Code %in% countries, "ISO"]
  }
  
  # Selected files for selected countries (1st/last time countries argument used)
  files.select <- lapply(files.all, function(y) sapply(countries, function(x) y[substr(y, 
  nchar(y)+config$input$cnt_part[1], nchar(y)+config$input$cnt_part[2])==tolower(x)])) 
  # no blanks, no home instrument, otherwise delete, see 4g function
  
  # Remove cases for no home instruments
  # only if home is specified
  if (!missing(home)) {
    if (!any(is.null(files.select[['ash']]))) {
      files.select[['ash']] <- files.select[['ash']][lapply(files.select[['ash']], length)!=0]
    }
    if (!any(is.null(files.select[['bsh']]))) {
      files.select[['bsh']] <- files.select[['bsh']][lapply(files.select[['bsh']], length)!=0]
    }
  }
  
  #
  # PREPARING
  #

  # Student achievement and background data, 
  # needed also if school is non-missing
  if (!missing(student) | !missing(school)) {
    if (!missing(student) & is.null(files.select[[config$input$student]])) {
      stop('cannot locate student data files')
    }
    
    suppressWarnings(suppressMessages(junk0 <- lapply(files.select[[config$input$student]], function(y) 
      read.spss(y, to.data.frame=TRUE, use.value.labels=use.labels))))
    student.data <- do.call('rbind', lapply(junk0, function(x) x[, c(config$input$student_colnames1,
      grep(config$input$student_pattern, names(x), value=TRUE), student, config$input$student_colnames2)]))
  }
  
  # Home background data
  if (!missing(home)) {
    if (is.null(files.select[[config$input$home]])) {
      stop('cannot locate home background data files')
    }
    
    suppressWarnings(suppressMessages(home.data <- do.call("rbind",                          # Merge [[2]] home
              lapply(files.select[[config$input$home]], function(y) 
                read.spss(y, to.data.frame=TRUE, use.value.labels=use.labels)[, c(           # Each dataset
                  config$input$home_colnames, home)]))))                                     # Selected
  }
  
  # School data
  if (!missing(school)) {
    if (is.null(files.select[[config$input$school]])) {
      stop('cannot locate school data files')
    }

    suppressWarnings(suppressMessages(school.data <- do.call("rbind",                      # Merge [[2]] school
           lapply(files.select[[config$input$school]], function(y) 
             read.spss(y, to.data.frame=T)[c(config$input$school_colnames, school)]))))    # Selected
  }
  
  # Teacher data
  if (!missing(teacher)) {
    if (is.null(files.select[[config$input$teacher[1]]]) | is.null(files.select[[config$input$teacher[2]]])) {
      stop('cannot locate teacher data files')
    }

    suppressWarnings(suppressMessages(teach.l <- do.call("rbind",                              
    lapply(files.select[[config$input$teacher[1]]], function(y) 
      read.spss(y, to.data.frame=TRUE, use.value.labels=use.labels)))))
    suppressWarnings(suppressMessages(teach.i <-  do.call("rbind",                              
    lapply(files.select[[config$input$teacher[2]]], function(y) 
      read.spss(y, to.data.frame=TRUE, use.value.labels=use.labels)[, c(
      config$input$teacher_colnames, teacher)]))))
  
    teacher.data <- merge(teach.l, teach.i, by=config$input$teacher_colnames)
  }

  #
  # MERGING
  #
  if (!missing(student)) {
    # student data is avaliable
    if (!missing(school)) {
      # + student + school
      if (!missing(teacher)) {
        # + student + school + teacher
        if (!missing(home)) {
          # + student + school + teacher + home
          student.data <- merge(student.data, home.data, by=config$input$student_ids, all.x=TRUE, suffixes=c(".st", ".hm"))
          student.data <- merge(teacher.data, student.data, by=config$input$student_ids)
          student.data <- student.data[, -c(grep("*.y", names(student.data)))]
          names(student.data) <- gsub("*.x", "", names(student.data))
          intsvy.all <- merge(student.data, school.data, by=config$input$school_ids, suffixes=c(".st", ".sc") )
          
        } else {
          # + student + school + teacher - home
          student.data <- merge(teacher.data, student.data, by=config$input$student_ids)
          student.data <- student.data[, -c(grep("*.y", names(student.data)))]
          names(student.data) <- gsub("*.x", "", names(student.data))
          intsvy.all <- merge(student.data, school.data, by=config$input$school_ids, suffixes=c(".st", ".sc"))
          
        }
      } else {
        # + student + school - teacher
        if (!missing(home)) {
          # + student + school - teacher + home
          student.data <- merge(student.data, home.data, by=config$input$student_ids, all.x=TRUE, suffixes=c(".st", ".hm"))
          intsvy.all <- merge(student.data, school.data, by=config$input$school_ids, suffixes=c(".st", ".sc"))
        } else {
          # + student + school - teacher - home
          intsvy.all <- merge(student.data, school.data, by=config$input$school_colnames[1:2], suffixes=c(".st", ".sc"))
        }
      }
    } else {
      # + student - school
      if (!missing(teacher)) {
        # + student - school + teacher
        if (!missing(home)) {
          # + student - school + teacher + home
          student.data <- merge(student.data, home.data, by=config$input$student_ids, all.x=TRUE, suffixes=c(".st", ".hm"))
          student.data <- merge(teacher.data, student.data, by=config$input$student_ids)
          student.data <- student.data[, -c(grep("*.y", names(student.data)))]
          names(student.data) <- gsub("*.x", "", names(student.data))
          intsvy.all <- student.data
          
        } else {
          # + student - school + teacher - home
          student.data <- merge(teacher.data, student.data, by=config$input$student_ids)
          student.data <- student.data[, -c(grep("*.y", names(student.data)))]
          names(student.data) <- gsub("*.x", "", names(student.data))
          intsvy.all <- student.data
        }
      } else {
        # + student - school - teacher
        if (!missing(home)) {
          # + student - school - teacher + home
          intsvy.all <- merge(student.data, home.data, all.x=TRUE, by=config$input$student_ids, suffixes=c(".st", ".hm"))
        } else {
          # + student - school - teacher - home
          intsvy.all <- student.data
        }
      }
    }
  } else {
    # student data is not avaliable
    if (!missing(school)) {
      # - student + school
      if (!missing(teacher)) {
        # - student + school + teacher
        if (!missing(home)) {
          # - student + school + teacher + home
          student.data <- merge(student.data, home.data, by=config$input$student_ids, suffixes=c(".st", ".hm"))
          student.data <- merge(teacher.data, student.data, by=config$input$student_ids)
          student.data <- student.data[, -c(grep("*.y", names(student.data)))]
          names(student.data) <- gsub("*.x", "", names(student.data))
          intsvy.all <- merge(student.data, school.data, by=config$input$school_ids, suffixes=c(".st", ".sc"))
          
        } else {
          # - student + school + teacher - home
          intsvy.all <- merge(teacher.data, school.data, by=config$input$school_ids, suffixes=c(".st", ".sc") )
        }

      } else {
        # - student + school - teacher
        if (!missing(home)) {
          # - student + school - teacher + home
          student.data <- merge(student.data, home.data, by=config$input$student_ids, suffixes=c(".st", ".hm"))
          intsvy.all <- merge(student.data, school.data, by=config$input$school_ids, suffixes=c(".st", ".sc"))
        } else {
          # - student + school - teacher - home
          intsvy.all <- merge(student.data, school.data, by=config$input$school_ids, suffixes=c(".st", ".sc") )
        }
      }
    } else {
      # - student - school
      if (!missing(teacher)) {
        # - student - school + teacher
        if (!missing(home)) {
          # - student - school + teacher + home
          student.data <- merge(student.data, home.data, by=config$input$student_ids, suffixes=c(".st", ".hm"))
          student.data <- merge(teacher.data, student.data, by=config$input$student_ids)
          student.data <- student.data[, -c(grep("*.y", names(student.data)))]
          names(student.data) <- gsub("*.x", "", names(student.data))
          intsvy.all <- student.data
        } else {
          # - student - school + teacher - home
          intsvy.all <- teacher.data
        }
        
      } else {
        # - student - school - teacher
        if (!missing(home)) {
          # - student - school - teacher + home
          intsvy.all <- merge(student.data, home.data, by=config$input$student_ids, suffixes=c(".st", ".hm"))
        } else {
          # - student - school - teacher - home
          # should not happen
        }
      }
    }
  }
  
  # Create country label variable (not possible to add labels to numeric factor, see to do list)
  intsvy.all$IDCNTRYL <- factor(intsvy.all$IDCNTRY,  
                               levels=country.list[country.list$Code %in% unique(intsvy.all$IDCNTRY), "Code"] ,        
                               labels= country.list[country.list$Code %in% unique(intsvy.all$IDCNTRY), "Country"])
  
  # table(student.data$IDCNTRYL, student.data$IDCNTRY) test of equality!

  return(droplevels(intsvy.all))
}
