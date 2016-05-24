################################################################################
# TODO LIST
# TODO: make 'exact' a vector.

################################################################################
# CHANGE LOG (last 20 changes)
# 09.01.2016: Added more attributes to result.
# 30.11.2015: Added attributes to result.
# 30.11.2015: Added parameter 'what' to specify columns to add.
# 15.12.2014: Changed parameter names to format: lower.case
# 30.09.2013: Fixed bug when exact=FALSE
# 17.09.2013: Updated example to support new 'getKit' structure.
# 31.07.2013: Added parameter 'ignore.case'.
# 28.05.2013: Fixed bug any(..., na.rm=TRUE)
# 21.05.2013: Added a second 'by' level and bugs fixed.
# 20.05.2013: Handle keys with no match.

#' @title Adds New Data Columns to a Data Frame
#'
#' @description
#' Adds values from columns in 'new.data' to 'data' by keys.
#'
#' @details
#' Information in columns in data frame 'new.data' is added to data frame
#' 'data' based on primary key value in column 'by.col', and optionally on
#' secondary key values in column 'then.by.col'.
#'   
#' @param data Data frame containing your main data.
#' @param new.data Data frame containing information you want to add to 'data'.
#' @param by.col character, primary key column.
#' @param then.by.col character, secondary key column.
#' @param exact logical, TRUE matches keys exact.
#' @param ignore.case logical, TRUE ignore case.
#' @param what character vector defining columns to add. Default is all new columns.
#' @param debug logical indicating printing debug information.
#' 
#' @return data.frame the original data frame containing additional columns.
#' 
#' @export
#' 
#' @examples
#' # Get marker names and alleles for Promega PowerPlex ESX 17.
#' x <- getKit("ESX17", what="Allele")
#' # Get marker names and colors for Promega PowerPlex ESX 17.
#' y <- getKit("ESX17", what="Color")
#' # Add color information to allele information.
#' z <- addData(data=x, new.data=y, by.col="Marker")
#' print(x)
#' print(y)
#' print(z)

addData <- function(data, new.data, by.col, then.by.col=NULL, exact=TRUE,
                    ignore.case=TRUE, what=NULL, debug=FALSE){
  
  # Parameters that are changed by the function must be saved first.
  attr_data <- substitute(data)

  if(debug){
    print(paste("IN:", match.call()[[1]]))
  }
  
  # Prepare -------------------------------------------------------------------
  
  # Remove unused colums.
  if(!is.null(what)){
    new.data <- new.data[ , c(by.col, then.by.col, what)]
  }

  # Get column names.
  colNames <- names(data)
  colNamesNew <- names(new.data)
  colNamesNew <- colNamesNew[colNamesNew!=by.col]
  if(!is.null(then.by.col)){
    colNamesNew <- colNamesNew[colNamesNew!=then.by.col]
  }
  colNamesNew <- colNamesNew[!colNamesNew %in% colNames]
  
  # Add new columns to data.
  for(c in seq(along=colNamesNew)){
    data[colNamesNew[c]] <- NA
  }
  
  # Get unique identifiers.
  keysData <- unique(as.character(data[ , by.col]))
  keysNew <- unique(as.character(new.data[ , by.col]))

  # Select keys.
  if(exact){
    keys <- keysData
  } else {
    keys <- keysNew
  }
  
  if(debug){
    print("keys")
    print(keys)
  }
  
  # Loop through keys.
  for(k in seq(along=keys)){
    
    # Select rows.
    if(exact){
      selectedData <- data[ , by.col] == keys[k]
      selectedNewData <- new.data[ , by.col] == keys[k]
    } else {
      selectedData <- grepl(keys[k], data[ , by.col],
                            ignore.case=ignore.case)
      selectedNewData <- grepl(keys[k], new.data[ , by.col],
                               ignore.case=ignore.case)
    }

    if(debug){
      print("selectedData")
      print(data[selectedData, ])
      print("selectedNewData")
      print(new.data[selectedNewData, ])
    }

    if(!is.null(selectedData)){
      
      if(is.null(then.by.col)){
        
        if(any(selectedData, na.rm = TRUE) && any(selectedNewData, na.rm = TRUE)){
          for(c in seq(along=colNamesNew)){
            
            dataLen <- length(data[selectedData , colNamesNew[c]])
            dataNewLen <- length(new.data[selectedNewData , colNamesNew[c]])
            
            if(dataLen == dataNewLen){
              
              # Add new data.
              data[selectedData , colNamesNew[c]] <- 
                new.data[selectedNewData , colNamesNew[c]]
            
            } else {
              
              uniqueNewData <- unique(new.data[selectedNewData , colNamesNew[c]])
              
              # See if all the same value.
              if(length(uniqueNewData) == 1) {
                # Add new data.
                data[selectedData , colNamesNew[c]] <- uniqueNewData
              } else {
                warning(paste("Ambiguous data could not be added at key",
                              keys[k]),
                        call. = TRUE)
              }
            }
          }
        }
        
      } else {
        
        # Get unique identifiers.
        keysData2 <- unique(as.character(data[ , then.by.col]))
        keysNew2 <- unique(as.character(new.data[ , then.by.col]))
        
        # Select keys.
        if(exact){
          keys2 <- keysData2
        } else {
          keys2 <- keysNew2
        }

        if(debug){
          print("keys2")
          print(keys2)
        }

        # Loop through keys.
        for(k2 in seq(along=keys2)){
          
          # Select rows.
          if(exact){
            selectedData2 <- data[ , then.by.col] == keys2[k2] & selectedData
            selectedNewData2 <- new.data[ , then.by.col] == keys2[k2] & selectedNewData
          } else {
            selectedData2 <- grepl(keys2[k2], data [ , then.by.col],
                                   ignore.case=ignore.case) & selectedData
            selectedNewData2 <- grepl(keys2[k2], new.data[ , then.by.col],
                                      ignore.case=ignore.case) & selectedNewData
          }
          
          if(debug){
            print("selectedData2")
            print(data[selectedData2, ])
            print("selectedNewData2")
            print(new.data[selectedNewData2, ])
          }

          if(any(selectedData2, na.rm = TRUE) && any(selectedNewData2, na.rm = TRUE)){
            
            for(c2 in seq(along=colNamesNew)){
              
              dataLen <- length(data[selectedData2 , colNamesNew[c2]])
              dataNewLen <- length(new.data[selectedNewData2 , colNamesNew[c2]])
              
              if(dataLen == dataNewLen){
                
                if(debug){
                  print("SelectedData")
                  print(data[selectedData2 , colNamesNew[c2]])
                  print("SelectedNewData")
                  print(new.data[selectedNewData2 , colNamesNew[c2]])
                }

                # Add new data.
                data[selectedData2 , colNamesNew[c2]] <- 
                new.data[selectedNewData2 , colNamesNew[c2]]
              
              } else {
                
                uniqueNewData <- unique(new.data[selectedNewData2 , colNamesNew[c2]])
                
                if(debug){
                  print("uniqueNewData")
                  print(uniqueNewData)
                }

                # See if all the same value.
                if(length(uniqueNewData) == 1) {
                  # Add new data.
                  data[selectedData2 , colNamesNew[c2]] <- uniqueNewData
                  
                } else {
                  warning(paste("Ambiguous data could not be added at key",
                                keys[k],"sub key:", keys2[k2]),
                       call. = TRUE)
                }
              }
            }
          }
        }
      }
    }
  }
  
  # Add attributes to result.
  attr(data, which="addData, strvalidator") <- as.character(utils::packageVersion("strvalidator"))
  attr(data, which="addData, call") <- match.call()
  attr(data, which="addData, date") <- date()
  attr(data, which="addData, data") <- attr_data
  attr(data, which="addData, new.data") <- substitute(new.data)
  attr(data, which="addData, by.col") <- by.col
  attr(data, which="addData, then.by.col") <- then.by.col
  attr(data, which="addData, exact") <- exact
  attr(data, which="addData, ignore.case") <- ignore.case
  attr(data, which="addData, what") <- what

  if(debug){
    print(paste("EXIT:", match.call()[[1]]))
  }
  
  return(data)
  
}
