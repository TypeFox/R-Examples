# 
# Parse the formula given to rollply
# 

parse_formula <- function(form, enclos) { 
  
  # If it has been already parsed, then we bail
  if ( ! is.null(attr(form, "parsed")) ) { 
    return(form)
  }
  
  # Textual form of input
  text <- gsub(' ', '', deparse(form))
  
  # Handle the case when group delimiter | is not specified
  if ( ! grepl("|",text, fixed=TRUE)) 
    text <- paste0(text,"|")
  
  # Discard stuff before ~ 
  # left_side  <- str_extract(text,perl("^.*(?=~)"))
  right_side <- stringr::str_extract(text, stringr::regex("(?<=~).*(?=\\|)"))
  groups     <- stringr::str_extract(text, stringr::regex("(?<=\\|).*$"))
  
  
  # Process text extracted from formula
  vars <- list(right_side, groups)
  vars <- lapply(vars, function(str) unlist(strsplit(str, "+", fixed=TRUE)))
  vars <- lapply(vars, function(str) str[!grepl("\\+",str)] )
  vars <- lapply(vars, as.quoted_or_NA, enclos)
  names(vars) <- c("coords","groups")
  
  attr(vars, "parsed") <- TRUE
  
  return(vars)
}

as.quoted_or_NA <- function(textvec, enclos) { 
  
  if (length(textvec) == 0) { 
    return(NA)
  } else { 
    return(plyr::as.quoted(textvec, env=enclos))
  }
}

