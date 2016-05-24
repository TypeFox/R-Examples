################################################################################
# TODO LIST
# TODO: ...

################################################################################
# CHANGE LOG (last 20 changes)
# 09.01.2016: Added attributes to result.
# 29.08.2015: Added importFrom.
# 25.05.2015: Corrected parameter description.
# 11.05.2015: Accepts (the first) column name containing the string 'Sample'
# as alternative to colum name 'Sample.Name'. 'Sample' is case in-sensitive.
# 04.05.2015: Added 'Sample.File.Name' as a defined alternative to Sample.Name.
# 15.12.2014: Changed parameter names to format: lower.case
# 28.04.2014: More robust and handles '+' and '-' in sample names.
# 14.01.2014: Support dataframes without a 'Sample.Name' column.
# 27.10.2013: Fixed bug when 'samples'=NULL and 'invert.s'=TRUE.
# 27.04.2013: Fixed error when result is only 1 column.
# <27.04.2013: Roxygenized.
# <27.04.2013: new name 'trimData' -> 'trim'
# <27.04.2013: remove NA/empty cols.
# <27.04.2013: handle no sample/ no columns.

#' @title Trim Data
#'
#' @description
#' Extract data from a dataset.
#'
#' @details
#' Simplifies extraction of specific data from a larger dataset.
#' Look for samples in column named 'Sample.Name', 'Sample.File.Name', or
#' the first column containing the string 'Sample' in mentioned order
#' (not case sensitive). Remove unwanted columns.
#' 
#' @param data data.frame with genotype data.
#' @param samples string giving sample names separated by pipe (|).
#' @param columns string giving column names separated by pipe (|).
#' @param word logical indicating if a word boundary should be added to 
#'  \code{samples} and \code{columns}.
#' @param ignore.case logical, TRUE ignore case in sample names.
#' @param invert.s logical, TRUE to remove matching samples from 'data',
#' FALSE to remove samples NOT matching (i.e. keep matching samples).
#' @param invert.c logical, TRUE to remove matching columns from 'data',
#' FALSE to remove columns NOT matching (i.e. keep matching columns).
#' while TRUE will remove columns NOT given.
#' @param rm.na.col logical, TRUE columns with only NA are removed from 'data'
#' while FALSE will preserve the columns.
#' @param rm.empty.col logical, TRUE columns with no values are removed from 'data'
#' while FALSE will preserve the columns.
#' @param missing value to replace missing values with.
#' @param debug logical indicating printing debug information.
#' 
#' @return data.frame with extracted result.
#' 
#' @export
#' 
#' @importFrom utils head
#' 


trim <- function(data, samples=NULL, columns=NULL, 
	word=FALSE, ignore.case=TRUE, invert.s=FALSE, invert.c=FALSE,
	rm.na.col=TRUE, rm.empty.col=TRUE, missing=NA, debug=FALSE){

  # Parameters that are changed by the function must be saved first.
  attr_data <- substitute(data)
  attr_columns <- columns
  
  # Variables.
  colNames <- columns
  
  if(debug){
    print(paste("IN:", match.call()[[1]]))
    print("data:")
    print(head(data))
    print("samples:")
    print(samples)
    print("columns:")
    print(columns)
    print("word:")
    print(word)
    print("ignore.case:")
    print(ignore.case)
    print("invert.s:")
    print(invert.s)
    print("invert.c:")
    print(invert.c)
    print("rm.na.col:")
    print(rm.na.col)
    print("rm.empty.col:")
    print(rm.empty.col)
    print("missing:")
    print(missing)
  }

  # Ignore case. NB! Must be before add word boundary.
  if(ignore.case){

    # Convert to upper case.
    samples <- toupper(samples)
    columns <- toupper(columns)
    
    if(length(samples) == 0){
      samples <- NULL
    }
    
    if(length(columns) == 0){
      columns <- NULL
    }
    
    if(debug){
      print("After ignore.case.")
      print("samples:")
      print(samples)
      print("columns:")
      print(columns)
    }
    
  }

  # Add word boundary. NB! Must be after ignore case.
  if(word){
    
    if(!is.null(samples)){
      samples <- gsub("|", "\\b|\\b", samples, fixed=TRUE)
      samples <- paste("\\b", samples, "\\b", sep="")
    }
    if(!is.null(columns)){
      columns <- gsub("|", "\\b|\\b", columns, fixed=TRUE)
      columns <- paste("\\b", columns, "\\b", sep="")
    }

    if(debug){
      print("After adding word boundary.")
      print("samples:")
      print(samples)
      print("columns:")
      print(columns)
    }
    
  }
  
  # Check for and escape '+' and '-' in sample names.
  if(any(grepl("+", samples, fixed=TRUE))){
    samples <- gsub("+", "\\+", samples, fixed=TRUE)
    message("'+' in sample names escaped")
  }
  if(any(grepl("-", samples, fixed=TRUE))){
    samples <- gsub("-", "\\-", samples, fixed=TRUE)
    message("'-' in sample names escaped")
  }
  
  # Grab samples --------------------------------------------------------------
  
  # Check if column 'Sample.Name' exist.
  if("Sample.Name" %in% names(data)){
    
    # Grab rows.
    if(ignore.case){
      sampleNames <- toupper(as.character(data$Sample.Name))
    } else {
      sampleNames <- as.character(data$Sample.Name)
    }

    if(debug){
      print("Pattern for samples")
      print(head(samples))
      print("String")
      print(head(sampleNames))
    }
    
    if(is.null(samples)){
      
      # Default is to keep all samples.
      rows <- rep(TRUE, length(sampleNames))
      
    } else {
      
      # Get matching rows.
      rows <- grepl(samples, sampleNames, fixed=FALSE)
      
      # Invert selection of samples.
      if(invert.s){
        rows <- !rows
      }
      
    }

  # Check if column 'Sample.File.Name' exist.
  } else if("Sample.File.Name" %in% names(data)){
    
    # Grab rows.
    if(ignore.case){
      sampleNames <- toupper(as.character(data$Sample.File.Name))
    } else {
      sampleNames <- as.character(data$Sample.File.Name)
    }
    
    if(debug){
      print("Pattern for samples")
      print(head(samples))
      print("String")
      print(head(sampleNames))
    }
    
    if(is.null(samples)){
      
      # Default is to keep all samples.
      rows <- rep(TRUE, length(sampleNames))
      
    } else {
      
      # Get matching rows.
      rows <- grepl(samples, sampleNames, fixed=FALSE)
      
      # Invert selection of samples.
      if(invert.s){
        rows <- !rows
      }
      
    }

  # Check if any column containing 'Sample' exist.
  } else if(any(grepl("SAMPLE", names(data), ignore.case=TRUE))){
    
    # Get (first) column name containing "Sample".
    sampleCol <- names(data)[grep("SAMPLE", names(data), ignore.case=TRUE)[1]]
    
      # Grab rows.
      if(ignore.case){
        sampleNames <- toupper(as.character(data[, sampleCol]))
      } else {
        sampleNames <- as.character(data[, sampleCol])
      }
      
      if(debug){
        print("Pattern for samples")
        print(head(samples))
        print("String")
        print(head(sampleNames))
      }
      
      if(is.null(samples)){
        
        # Default is to keep all samples.
        rows <- rep(TRUE, length(sampleNames))
        
      } else {
        
        # Get matching rows.
        rows <- grepl(samples, sampleNames, fixed=FALSE)
        
        # Invert selection of samples.
        if(invert.s){
          rows <- !rows
        }
        
      }
      
  } else {
    
    # Keep all rows.
    rows <- rep(TRUE, nrow(data))
    
  }
  
  if(debug){
    print(paste("Grab samples:", paste(sampleNames[rows], collapse=",")))
  }
  
  # Grab columns --------------------------------------------------------------
  
  if(ignore.case){
    columnNames <- toupper(names(data))
  } else {
    columnNames <- names(data)
  }

  if(debug){
    print("Pattern for columns")
    print(head(columns))
    print("String")
    print(head(columnNames))
  }
  
  
	if(is.null(columns)){
    
		# Default is to keep all columns.
		columns <- rep(TRUE, length(columnNames))
    
	} else {
    
    # Get matching columns.
		columns <- grepl(columns, columnNames, fixed=FALSE)
    
		# Invert selection of columns.
		if(invert.c){
		  columns <- !columns
		}
		
	}

  if(debug){
    print(columns)
    print(paste("Grab columns:", paste(names(data)[columns], collapse=",")))
  }
  
	# Trim data.
	data <- data[rows, columns]

	if(!is.null(missing)){
		data[data==""] <- missing
	}

  if(rm.empty.col){
    if(!is.null(ncol(data))){
		  data <- data[ , colSums(data=="") != nrow(data) | colSums(is.na(data))>0]
    }
	}

	if(rm.na.col){
	  if(!is.null(ncol(data))){
	    data <- data[,colSums(is.na(data))<nrow(data)]
	  }
	}

  if(debug){
    print(paste("EXIT:", match.call()[[1]]))
  }

  if(is.null(ncol(data))){
    data <- data.frame(Name=data)
    names(data) <- colNames
  }
	
	# Add attributes to result.
	attr(data, which="trim, strvalidator") <- as.character(utils::packageVersion("strvalidator"))
	attr(data, which="trim, call") <- match.call()
	attr(data, which="trim, date") <- date()
	attr(data, which="trim, data") <- attr_data
	attr(data, which="trim, samples") <- samples
	attr(data, which="trim, columns") <- attr_columns
	attr(data, which="trim, word") <- word
	attr(data, which="trim, ignore.case") <- ignore.case
	attr(data, which="trim, invert.s") <- invert.s
	attr(data, which="trim, invert.c") <- invert.c
	attr(data, which="trim, rm.na.col") <- rm.na.col
	attr(data, which="trim, rm.empty.col") <- rm.empty.col
	attr(data, which="trim, missing") <- missing
	
  return(data)
	
}
