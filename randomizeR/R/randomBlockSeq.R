#' @include rpbrSeq.R
#' @include rtbdSeq.R
NULL


###############################################
# --------------------------------------------#
# Class randomBlockSeq                              #
# --------------------------------------------#
###############################################

# Random Block Design class
# 
# @name randomBlockSeq
setClassUnion("randomBlockSeq", members = c("rRtbdSeq", "rRpbrSeq"))

# --------------------------------------------
# Show function for randSeq
# --------------------------------------------

setMethod("show", "randomBlockSeq", function(object) {
  # headline
  cat("\nObject of class \"", class(object)[1],"\"\n\n", sep = "")
  # crop the method from the class name of the randPar object
  cat("design =", getDesign(object), " \n")  
  # iterate through all slots of the object
  names <- slotNames(object)
  names <- names[!(names %in% c("M", "bc"))] # without M and bc
  for(name in names) {
    cat(name, "=", slot(object, name), "\n")
  }  
  
  if (dim(object@M)[1] <= 3) {
    # The matrix M is printed seperately dependent on its size.
    # sequences <- getRandList(object)
    if (object@N <= 8) {
      sequences <- apply(getRandList(object), 1, function(x) paste(x, collapse = " "))
      # bc <- blocks(object)
      blockConst <- unlist(lapply(blocks(object), function(x) paste(x, collapse = " ")))
    } else {
      sequences <- apply(getRandList(object), 1, function(x)
        paste(c(x[1:8], "..."), collapse=" "))
      
      blockConst <- unlist(lapply(blocks(object), function(x) {
        if (length(x) <= 5) {  
          return(paste(x, collapse = " "))
        } else {
          return(paste(c(x[1:5], "..."), collapse = " "))
        }
      }))
    }
    
    print.data.frame <- function(m) {
      write.table(format(m, justify = "left"),
                  row.names = F, col.names = T, quote = F)
    }
    cat("\n") 
    print(data.frame("RandomizationSeqs" = paste(sequences, "\t"),
                     BlockConst = blockConst))
  } else {
    # The matrix M is printed seperately dependent on its size.
    # sequences <- getRandList(object)
    numberSeq <- dim(object@M)[1] 
    object@M <- object@M[1:4, ]  
    if (object@N <= 8) {
      sequences <- apply(getRandList(object)[1:3 , ], 1, function(x) paste(x, collapse = " "))
      # bc <- blocks(object)
      blockConst <- unlist(lapply(object@bc[1:3], function(x) paste(x, collapse = " ")))
    } else {
      sequences <- apply(getRandList(object)[1:3 , ], 1, function(x)
        paste(c(x[1:8], "..."), collapse = " "))
      
      
      blockConst <- unlist(lapply(object@bc[1:3], function(x) {
        if (length(x) <= 5) {  
          return(paste(x, collapse = " "))
        } else {
          return(paste(c(x[1:5], "..."), collapse = " "))
        }
      }))
    }
    
    sequences <- c(sequences, "...")
    blockConst <- c(blockConst, "...")
    cat("\nThe first 3 of", numberSeq, "sequences of M:\n")
    
    print.data.frame <- function(m) {
      write.table(format(m, justify = "left"),
                  row.names = F, col.names = T, quote = F)
    }
    
    cat("\n") 
    print(data.frame("RandomizationSeqs" = paste(sequences, "\t"),
                     BlockConst = blockConst))
    
  }  
})