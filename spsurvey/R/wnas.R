wnas <- function(data) {

################################################################################
# Function: wnas
# Programmer: Tom Kincaid
# Date: November 16, 2000
# Last Revised: May 5, 2006
# Description:
#   This function removes missing values from data for which the mode is one of 
#   the following: numeric, logical, or character.  Data that is not one of 
#   those modes will cause the function to terminate with an error message.  
#   For numeric data this function removes values that are not finite, i.e., 
#   missing value (NA), not a number (NaN), infinity (Inf), and minus infinity
#   (-Inf).  For logical data this function removes missing values (NA).  For 
#   character data the following values are removed: "", "NA", NA (R only), 
#   "NaN", "Inf", and "-Inf".  For a factor this function removes the following 
#   values: NA, NaN, Inf, and -Inf.  For a vector this function returns the 
#   vector with the indicated values removed.  For a data frame this function 
#   returns the data frame with rows removed that contain at least one indicated 
#   value.  For a list with components that are the same length, the list is 
#   converted to a data frame.  For a list with components that are not the same 
#   length, the function prints an error message and terminates.  When the 
#   process of removing missing values produces an object that no longer
#   contains any elements (vector) or rows (data frame or list), then a NULL
#   object is returned.
################################################################################

   if(is.list(data)) {
      if(any(match(sapply(data, mode), c("numeric", "logical", "character"), nomatch=0) == 0))
         stop("\nAt least one component in the list passed to the missing value function was not \none of the following modes: numeric, logical, or character.")
      if(length(unique(sapply(data, length))) > 1)
         stop("\nComponents in the list passed to the missing value function must be the same \nlength.")
      if(!is.data.frame(data))
         data <- as.data.frame(data)
      n <- length(data)
      for (i in 1:n) {
         if(all(is.na(data[,i]))) {
            data <- NULL
            return(data)
         } else if(is.numeric(data[,i]))
            data <- data[match(data[,i], c(NA, NaN, Inf, -Inf), nomatch=0) == 0,]
         else if(is.logical(data[,i]))
            data <- data[!is.na(data[,i]),]
         else if(is.character(data[,i]))
            data <- data[match(data[,i], c("", "NA", NA, "NaN", "Inf", "-Inf"), nomatch=0) == 0,]
         else if(is.factor(data[,i]))
            data <- data[match(data[,i], c(NA, NaN, Inf, -Inf), nomatch=0) == 0,]
         else
            stop("\nThe input data was not recognized by the missing value function.")
      }
      if(nrow(data) > 0) {
         for(i in 1:n) {
            if(is.factor(data[,i]))
               data[,i] <- factor(data[,i])
         }
      }
   } else {
      if(match(mode(data), c("numeric", "logical", "character"), nomatch=0) == 0)
         stop("\nData passed to the missing value function was not one of the following modes: \nnumeric, logical, or character.")
      else if(all(is.na(data)))
         data <- NULL
      else if(is.numeric(data))
            data <- data[match(data, c(NA, NaN, Inf, -Inf), nomatch=0) == 0]
      else if(is.logical(data))
         data <- data[!is.na(data)]
      else if(is.character(data))
         data <- data[match(data, c("", "NA", NA, "NaN", "Inf", "-Inf"), nomatch=0) == 0]
      else if(is.factor(data))
         data <- factor(data[match(data, c(NA, NaN, Inf, -Inf), nomatch=0) == 0])
      else
         stop("\nThe input data was not recognized by the missing value function.")
   }
   data
}

