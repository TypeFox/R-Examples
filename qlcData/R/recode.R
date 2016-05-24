# ======================
# recode data according to specifications in recoding
#=======================

recode <- function(data,recoding) {

  # expand the possible shortcuts in the formulation of a recoding
  recoding <- read.recoding(recoding)

  # recoding of a single new attribute
  makeAttribute <- function(recoding) {

    # when doNotRecode is specified, do not recode attributes
    if (!is.null(recoding$doNotRecode)) {
      newAttribute <- data[,recoding$doNotRecode, drop = FALSE]
    } else {
      
      recoding$link[recoding$link == 0]  <- NA
      
      # simple when it is based on a single old attribute
      if (length(recoding$recodingOf) == 1) {
      newAttribute <- data[,recoding$recodingOf, drop = FALSE]
      levels(newAttribute[,1]) <- recoding$values[recoding$link]
      colnames(newAttribute) <- recoding$attribute
      return(newAttribute)
      } else {
        
        # a bit more complex for combinations of attributes
        # this can probably be made more efficient!
        newAttribute <- data[,recoding$recodingOf, drop = FALSE]
        newAttribute <- apply(newAttribute,1,function(x){paste(x, collapse = " + ")})
        match <- expand.grid(
          sapply(recoding$recodingOf, function(x){ 
            c(levels(data[,x]),NA) 
            }, simplify = FALSE )
          )
        match <- apply(match,1,function(x){paste(x, collapse = " + ")})
        newAttribute <- factor(newAttribute, levels = match)
        levels(newAttribute) <- recoding$values[recoding$link]
        newAttribute <- as.data.frame(newAttribute)
        colnames(newAttribute) <- recoding$attribute
        return(newAttribute)
      }
    }
  }
  
  # Make the recoding and return result
  result <- as.data.frame(sapply(recoding, makeAttribute,simplify=F))
  return(result)
}

