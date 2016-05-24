##################################################################################
##                      METHODS FOR S3 GENERICS                                 ##
##################################################################################

## ---------------------------------------------------------------------------##
## DATA FRAME COERCION METHOD

## This is a method for the as.data.frame S3 generic. We need this to intercept the
## DRAC list object after it hast passed the actual list-method. After it was 
## coerced to a data.frame we assign new column names (DRAC ID keys) and 
## make sure that all columns are either of class 'character' or 'numeric'.
## Finally, we attach a further class name to identify it as a valid DRAC object 
## when passed to use_DRAC

#' @export
as.data.frame.DRAC.list <- function(x, row.names = NULL, optional = FALSE, ...) {
  DF <- as.data.frame.list(x)
  colnames(DF) <- paste0("TI:", 1:ncol(DF))
  for (i in 1:ncol(DF)) {
    if (is.factor(DF[ ,i])) 
      DF[ ,i] <- as.character(DF[, i])
  }
  class(DF) <- c("data.frame", "DRAC.data.frame")
  return(DF)
}


## ---------------------------------------------------------------------------##
## PRINT METHOD

#' @export
print.DRAC.highlights <- function(x, ...) {
  x <- as.list(x)
  names <- names(x)
  mapply(function(el, name) { 
    cat(paste0(attributes(el)$key, " = ", name,":\n  ", paste(el, collapse = ",\n  "), "\n"))
    }, x, names)
}

#' @export
print.DRAC.list <- function(x, blueprint = FALSE, ...) {
  
  ## CASE 1: Pretty print the structure of the DRAC list
  if (!blueprint) {
    limit <- 80
    
    for (i in 1:length(x)) {
      # for pretty printing we insert newlines and tabs at specified lengths
      ls <- attributes(x[[i]])$description
      ls.n <- nchar(ls)
      ls.block <- floor(ls.n / limit)
      strStarts <- seq(0, ls.n, limit)
      strEnds <- seq(limit-1, ls.n + limit, limit)
      blockString <- paste(mapply(function(start, end) { 
        trimmedString <- paste(substr(ls, start, end), "\n\t\t\t")
        if (substr(trimmedString, 1, 1) == " ")
          trimmedString <- gsub("^[ ]*", "", trimmedString)
        return(trimmedString)
      }, strStarts, strEnds), collapse="")
      
      msg <- paste(attributes(x[[i]])$key, "=>",names(x)[i], "\n",
                   "\t VALUES =", paste(x[[i]], collapse = ", "), "\n",
                   "\t ALLOWS 'X' = ", attributes(x[[i]])$allowsX, "\n",
                   "\t REQUIRED =", attributes(x[[i]])$required, "\n",
                   "\t DESCRIPTION = ", blockString, "\n"
      )
      if (!is.null(levels(x[[i]]))) {
        msg <- paste(msg,
                     "\t OPTIONS = ", paste(levels(x[[i]]), collapse = ", "),
                     "\n\n")
      } else {
        msg <- paste(msg, "\n")
      }
      cat(msg)
    }
  }
  
  ## CASE 2: Return a 'blueprint' that can be copied from the console to a
  ## script so the user does not need to write down all >50 fields by hand
  if (blueprint) {
    var <- as.list(sys.call())[[2]]
    names <- names(x)
    
    for (i in 1:length(x)) {
      
      # in case of factors also show available levels as comments so you don't
      # have to look it up
      if (is.factor(x[[i]]))
        options <- paste("# OPTIONS:", paste(levels(x[[i]]), collapse = ", "))
      else
        options <- ""
      
      # determine if values need brackets (strings)
      if (is.numeric(x[[i]]) | is.integer(x[[i]]))
        values <- paste(x[[i]], collapse = ", ")
      if (is.character(x[[i]]) | is.factor(x[[i]]))
        values <- paste0("'", paste0(x[[i]], collapse = "', '"), "'")
      
      cat(paste0(var, "$`", names[i], "` <- c(", values,") ", options ,"\n"))
    }
    message("\n\t You can copy all lines above to your script and fill in the data.")
  }
}


## ---------------------------------------------------------------------------##
## DOUBLE SQUARE BRACKETS METHOD

#' @export
`[[<-.DRAC.list` <- function(x, i, value) {
  
  
  ## REJECT ALL INADEQUATE CLASSES ----
  acceptedClasses <- c("integer", "character", "numeric", "factor")
  if (is.na(match(class(value), acceptedClasses))) {
    warning(paste("I cannot use objects of class", class(value)), 
            call. = FALSE)
    return(x)
  }
  
  ## CHECK INPUT LENGTH ----
  length.old <- length(x[[i]])
  length.new <- length(value)
  
  if (length.old != length.new) {
    warning(paste(names(x)[i], ": Input must be of length", length.old), 
            call. = FALSE)
    return(x)
  }
  
  ## CHECK INPUT CLASS ----
  class.old <- class(x[[i]])
  class.new <- class(value)
  
  ## CHECK INPUT FIELDS THAT ALLOW 'X' -----
  # the following checks apply to fields that are normally numeric, but also 
  # accept 'X' as input. this EXCLUDES factors!
  if (class.old != "factor") {
    # some input fields allow 'X' as input, so in terms of R can be of class
    # "character" or "numeric/integer". hence, we check if input is "X" and 
    # if the filed allows it. If so, we change the old class to "character".
    if (any(value == "X") && attributes(x[[i]])$allowsX) {
      
      if (any(is.na(as.numeric(value[which(value != "X")])))) {
        warning(paste("Cannot coerce <", value[which(value != "X")], "> to a numeric value.",
                      "Input must be numeric or 'X'."), 
                call. = FALSE)
        return(x)
      }
      class.old <- "character" 
    }
    
    # where the input field is alreay "X" we have to check whether the new
    # non-character input is allowed
    if (any(x[[i]] == "X") && attributes(x[[i]])$allowsX) {
      if (any(is.na(as.numeric(value[which(value != "X")])))) {
        warning(paste("Cannot coerce <", value[which(value != "X")], "> to a numeric value.",
                      "Input must be numeric or 'X'. \n"), 
                call. = FALSE)
        return(x)
      }
      class.new <- "character"
      value <- as.character(value)
    }
    
    # when a numeric input field was inserted an "X" it was coerced to class
    # character. since we are now allowed to insert any character (later tests)
    # we need to make sure that the new input can be coerced to class numeric.
    # and if the new values are numeric, we coerce them to character
    if (attributes(x[[i]])$allowsX && class.old == "character") {
      if (any(is.na(as.numeric(value[which(value != "X")])))) {
        warning(paste("Cannot coerce <", value[which(value != "X")], "> to a numeric value.",
                      "Input must be numeric or 'X'. \n"), 
                call. = FALSE)
        return(x)
      } 
      class.new <- "character"
      value <- as.character(value)
    }
  }
  
  
  # numeric input can be both of class 'integer' or 'numeric'. We will
  # allow any combination and reject only non-numeric/integer input
  if (class.old == "numeric" || class.old == "integer") {
    if (class.new != "numeric" && class.new != "integer") {
      warning(paste(names(x)[i], ": Input must be of class", class.old),
              call. = FALSE)
      return(x)
    }
  }
  
  # for 'factor' and 'character' elements only 'character' input is allowed 
  if (class.old == "factor" || class.old == "character") {
    if (class.new != "character") {
      warning(paste(names(x)[i], ": Input must be of class", "character"),
              call. = FALSE)
      return(x)
    }
  }
  
  ## CHECK IF VALID OPTION ----
  # in case of 'factor's the user is only allowed input that matches one of 
  # the options specified by the factor levels. if it is a valid option,
  # the input is converted to a factor to keep the information.
  if (class.old == "factor") {
    levels <- levels(x[[i]])
    if (any(`%in%`(value, levels) == FALSE)) {
      warning(paste(names(x)[i], ": Invalid option. Valid options are:", paste(levels, collapse = ", ")),
              call. = FALSE)
      return(x)
    } else {
      value <- factor(value, levels)
    }
  }
  
  ## WRITE NEW VALUES ----
  # we strip our custom class and the attributes, pass the object to the default generic and 
  # finally re-attach our class and attributes
  tmp.attributes <- attributes(x[[i]])[names(attributes(x[[i]])) != "class"]
  class(x) <- "list"
  x <- `[[<-`(x, i, value)
  attributes(x[[i]]) <- tmp.attributes
  if (class.old == "factor")
    class(x[[i]]) <- "factor"
  class(x) <- c("DRAC.list", "list")
  return(x)
}

## ---------------------------------------------------------------------------##
## SINGLE SQUARE BRACKET METHOD

#' @export
`[<-.DRAC.list` <- function(x, i, value) {
  return(`[[<-`(x, i, value))
}

## ---------------------------------------------------------------------------##
## DOLLAR SIGN METHOD

#' @export
`$<-.DRAC.list`<- function(x, name, value) {
  # this is straightforward; retrieve the index and pass the object
  # to the custom [[<- function, which does the data verification
  index <- which(names(x) == name)
  x[[index]] <- value
  return(x)
}