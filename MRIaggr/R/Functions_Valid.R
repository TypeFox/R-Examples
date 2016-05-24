#### validation functions ####

validCharacter <- function(value, name = as.character(substitute(value)), validLength, 
                           validValues = "character", refuse.NULL = TRUE, method){
  
  if(is.null(value)){
    
    if(refuse.NULL == TRUE){
      stop(method, " : wrong specification of \'", name, "\' \n", 
           "\'", name, "\' must not be NULL \n")
    }
    
  }else{
    
    #### check size
    n.value <- length(value)
    
    if(!is.null(validLength) && n.value %in% validLength == FALSE){
      stop(method, " : wrong specification of \'", name, "\' \n", 
           "\'", name, "\' must have length ", paste(validLength, collapse = " or "), "  \n", 
           "length(", name, ") : ", n.value, "\n")
    }
    
    #### check values
    
    if(identical(validValues,"character")){
      
      if(any(is.character(value) == FALSE)){
        stop(method, " : wrong specification of \'", name, "\' \n", 
             "\'", name, "\' must be a ", if(n.value == 1){"character"}else{"vector of characters"}," \n", 
             "is(", name, ") : ", paste(is(value), collapse = " "), "\n")
      }
      
    } else if(identical(validValues,"character_or_logical")){
      
      if(any( (is.character(value) == FALSE) * (is.logical(value) == FALSE) > 0 )){
        stop(method, " : wrong specification of \'", name, "\' \n", 
             "\'", name, "\' must be a ", if(n.value == 1){"character or logical"}else{"vector of characters or logicals"}," \n", 
             "is(", name, ") : ", paste(is(value), collapse = " "), "\n")
      }
      
    } else if(!is.null(validValues) && any(value %in% validValues == FALSE)){
      
      stop(method, " : wrong specification of \'", name, "\' \n", 
           "valid values for \'", name, "\' : ", if(refuse.NULL == FALSE){"NULL"}, " \"", paste(validValues, collapse = "\" \""), "\" \n", 
           "refused value",if(sum(value %in% validValues == FALSE)>1){"s"}," for \'", name, "\' : \"", paste(value[value %in% validValues == FALSE], collapse = " "), "\"\n")
      
    }
    
  }
  
}

validClass <- function(value, name = as.character(substitute(value)), validClass, 
                       superClasses = TRUE, method){
  
  if(superClasses == TRUE){
    
    if( all(is(value) %in% validClass == FALSE) ){
      stop(method, " : wrong specification of \'", name, "\' \n", 
           "superclass of \'", name, "\' must be one of the following \"", paste(validClass,collapse="\" \""), "\"  \n", 
           "proposed superclass : \"", paste(is(value),sep="\" \""), "\" \n")
    }  
    
  }else{
 
    if( class(value) %in% validClass == FALSE){
      stop(method, " : wrong specification of \'", name, "\' \n", 
           "class of \'", name, "\' must be \"", paste(validClass,collapse="\" \""),"\"  \n", 
           "proposed class : ", class(value)[[1]], "\n")
    }  
    
  }
  
  
  
}

validDim_vector <- function(value1, value2, name1 = as.character(substitute(value1)), name2 = as.character(substitute(value2)),
                            type = "length", method){
  
  #### initialisation
  # value1
  if (is.list(value1)) {value1 <- unlist(value1)} 
  else if (is.vector(value1)) {value1 <- length(value1)}
  
  # value2
  if (is.list(value2)) {value2 <- unlist(value2)}
  else if (is.vector(value2)) {value2 <- length(value2)}
  else if (is.matrix(value2)) {
    if (type == "nrow") {value2 <- nrow(value2)}
    if (type == "ncol") {value2 <- ncol(value2)}
  }

  #### tests
  if(value1 != value2){
    
    if(!is.null(name2)){
      stop(method, " : dimension mismatch between argument \'", name1, "\' and argument \'", name2, "\' \n", 
           "length(", name1, ") = ", value1[1], " \n", 
           type,"(", name2, ") = ", value2[1], " \n")
      
    }else{
      
      stop(method, " : wrong specification of argument \'", name1,"\' \n", 
           "\'", name1,"\' must have size ",value2[1]," \n", 
           "length(", name1, ") = ", value1[1], " \n")  
      
    }
    
  }
  
}

validDim_matrix <- function(value1, value2, name1 = as.character(substitute(value1)), name2 = as.character(substitute(value2)),
                            type = "both", method){

  if ( any(is(value1) %in% c("matrix", "Matrix", "dgCMatrix", "data.frame")) ) {
    value1 <- dim(value1)
  }
  
  if ( any(is(value2) %in% c("matrix", "Matrix", "dgCMatrix", "data.frame")) ) {
    value2 <- dim(value2)
  } else if(length(value2) == 1 && type == "ncol"){
    value2 <- c(NA,value2)
  }

  if(type %in% c("both","nrow") && value1[1] != value2[1]){
    
    if(!is.null(name2)){
      stop(method, " : dimension mismatch between argument \'", name1, "\' and argument \'", name2, "\' \n", 
           "nrow(", name1, ") = ", value1[1], " \n", 
           "nrow(", name2, ") = ", value2[1], " \n")
    }else{
      stop(method, " : wrong specification of argument \'", name1,"\' \n", 
           "\'", name1,"\' must have ",value2[1]," rows \n", 
           "nrow(", name1, ") = ", value1[1], " \n")  
    }
    
  }
  
  if(type %in% c("both","ncol") && value1[2] != value2[2]){
    
    if(!is.null(name2)){
      stop(method, " : dimension mismatch between argument \'", name1, "\' and argument \'", name2, "\' \n", 
           "ncol(", name1, ") = ", value1[2], " \n", 
           "ncol(", name2, ") = ", value2[2], " \n")
    }else{
      stop(method, " : wrong specification of argument \'", name1,"\' \n", 
           "\'", name1,"\' must have ",value2[2]," columns \n", 
           "ncol(", name1, ") = ", value1[2], " \n")  
    }
    
  }
  
}

validInteger <- function(value, name = as.character(substitute(value)), validLength, 
                         validValues = NULL, min = NULL, max = NULL, 
                         refuse.NA = TRUE, refuse.NULL = TRUE, refuse.duplicates = FALSE, method){
  
  validNumeric(value = value, name = name, validLength = validLength, min = min, max = max, 
               refuse.NA = refuse.NA, refuse.NULL = refuse.NULL, refuse.duplicates = refuse.duplicates, method = method)
  
  #### check integer
  if(any(value %% 1 > 0)){
    stop(method, " : wrong specification of \'", name, "\' \n", 
         "\'", name, "\' must contain integers not doubles \n",        
         "invalid value(s) in ", name, " : ", paste(value[value %% 1 > 0], collapse = " "), "\n")
  }
  
}

validLogical <- function(value, name = as.character(substitute(value)), validLength, 
                         refuse.NULL = TRUE, refuse.NA = TRUE, method){
  
  
  if(is.null(value)){
    
    #### NULL
    if(refuse.NULL == TRUE){
      stop(method, " : wrong specification of \'", name, "\' \n", 
           "\'", name, "\' must be logical ",if(refuse.NA == FALSE){"or NA"}," and not NULL \n")
    }
    
  }else{ 
    
    #### Size
    if(!is.null(validLength) && length(value) %in% validLength == FALSE){
      stop(method, " : wrong specification of \'", name, "\' \n", 
           "\'", name, "\' must have length ", paste(validLength, collapse = " or "), "  \n", 
           "length(", name, ") : ", length(value), "\n")
    } 
    
    #### Type
    if(any(is.logical(value) == FALSE)){
      stop(method, " : wrong specification of \'", name, "\' \n", 
           "\'", name, "\' must be ", if(refuse.NULL == FALSE){"NULL or "}, if(refuse.NA == FALSE){"NA or "},"TRUE or FALSE \n",        
           "is(", name, ") : ", paste(is(value), collapse = " "), "\n")
    }
    
    if(refuse.NA == TRUE && any(is.na(value)) ){
      stop(method, " : wrong specification of \'", name, "\' \n", 
           "\'", name, "\' must be logical ",if(refuse.NULL == FALSE){"or NULL"}," and not NA \n")
    }
    
  }
  
}

validNames <- function(value, name = as.character(substitute(value)),
                       validLength = NULL, validValues = NULL, method){
  
  ## type
  if(is.matrix(value)){
    value <- colnames(value)
  }
  
  if(is.data.frame(value) || is.list(value)){
    value <- names(value)
  }
  
  ## tests
  if(is.null(value)){
    
    stop(method, " : wrong specification of \'", name, "\' \n", 
         "\'", name, "\' must not be NULL \n")
    
  }else{
    
    #### check size
    n.value <- length(value)
    
    if(!is.null(validLength) && n.value %in% validLength == FALSE){
      stop(method, " : wrong specification of \'", name, "\' \n", 
           "\'", name, "\' must have ", paste(validLength, collapse = " or ")," names  \n", 
           "length(names(", name, ")) : ", n.value, "\n")
    }
    
    #### check content
    
    if(!is.null(validValues) && any(value %in% validValues == FALSE)){
      
      stop(method, " : wrong specification of \'", name, "\' \n", 
           "valid names for \'", name, "\' : \"",paste(validValues, collapse = "\" \""),"\" \n", 
           "refused names : \"", paste(value[value %in% validValues == FALSE], collapse = " "), "\"\n")  
      
    }
    
    if(any(duplicated(value))){
      stop(method, " : wrong specification of \'", name, "\' \n", 
           name, " must not contain duplicated names \n", 
           "duplicated names : \"", paste(value[duplicated(value)], collapse = " "), "\"\n")  
    }
    
  }
  
}

validNumeric <- function(value, name = as.character(substitute(value)), validLength,
                         validValues = NULL , min = NULL, max = NULL,
                         refuse.NA = TRUE, refuse.NULL = TRUE, refuse.duplicates = FALSE, method){
  
  if(is.null(value)){
    
    if(refuse.NULL == TRUE){
      stop(method, " : wrong specification of \'", name, "\' \n", 
           "\'", name, "\' must not be NULL \n")
    }
    
  }else{
    
    #### check length
    if(!is.null(validLength) && length(value) %in% validLength == FALSE){
      stop(method, " : wrong specification of \'", name, "\' \n", 
           "\'", name, "\' must have length ", paste(validLength, collapse = " or "), "  \n", 
           "length(", name, ") : ", length(value), "\n")
    }
    
    #### check NA
    if(refuse.NA == TRUE && any(is.na(value))){
      stop(method, " : wrong specification of \'", name, "\' \n", 
           "\'", name, "\' must not contain NA \n", 
           "index of NA values : ", which(paste(is.na(value), collapse = " ")), "\n")
    }
    
    #### check numeric
    if(any( (is.numeric(value) == FALSE) * (is.na(value) == FALSE) )){
      stop(method, " : wrong specification of \'", name, "\' \n", 
           "\'", name, "\' must be a numeric \n",        
           "is(", name, ") : ", paste(is(value), collapse = " "), "\n")
    }
    
    #### check duplicates
    if(refuse.duplicates == TRUE && any(duplicated(value))){
      stop(method, " : wrong specification of \'", name, "\' \n", 
           "\'", name, "\' contains duplicated values : \n",        
           paste(value[duplicated(value)], collapse = " "), "\n")
    }
    
    #### check min value
    if(!is.null(min) && any(stats::na.omit(value) < min)){
      stop(method, " : wrong specification of \'", name, "\' \n", 
           "\'", name, "\' must be bigger than ", min, " \n",        
           "invalid value(s) in ", name, " : ", paste(value[stats::na.omit(value) < min], collapse = " "), "\n")
    }
    
    #### check max value
    if(!is.null(max) && any(stats::na.omit(value) > max)){
      stop(method, " : wrong specification of \'", name, "\' \n", 
           "\'", name, "\' must be smaller than ", max, " \n",        
           "invalid value(s) in ", name, " : ", paste(value[stats::na.omit(value) > max], collapse = " "), "\n")
    }
    
    #### check valid values
    if(!is.null(validValues) && any(value %in% validValues == FALSE)){
      
      stop(method, " : wrong specification of \'", name, "\' \n", 
           "valid values for \'", name, "\' : ", if(refuse.NULL == FALSE){"NULL"}, " \"", paste(validValues, collapse = "\" \""), "\" \n", 
           "refused value",if(sum(value %in% validValues == FALSE)>1){"s"}," for \'", name, "\' : \"", paste(value[value %in% validValues == FALSE], collapse = " "), "\"\n")
      
    }
  }
}

validPath <- function(value, name = as.character(substitute(value)),
                      method){

  if(!is.null(value)){
    try_path <- dir.exists(value)
    
    if(any(is(try_path) == "try-error")){
      stop(method, " : wrong specification of \'", name, "\' \n", 
           "proposed ", name, " : \"", value, "\" \n", 
           "error : ", paste(try_path, collapse = " "), 
           "current ", name, " : ", getwd(), "\n")
    }
    
    
    if(substr(value, start = nchar(value), stop = nchar(value)) != "/"){
      warning(method, " : possible bad specification of \'", name, "\' \n", 
              "\'", name, "\' should end with a fsep (e.g. \"/\") \n", 
              "proposed ", name, " : ", value, "\n")
    }
  }
}
