##
## show-method for objects of class NeosAns
##
setMethod("show", "NeosAns", function(object){
  if(is.character(object@ans) && length(object@ans) == 1){
    cat("\n")
    cat(object@ans)
    cat("\n")
  } else {
    print(object@ans)
  }
})
##
## show-method for objects of class NeosXml
##
setMethod("show", "NeosXml", function(object){
  print(object@xml)
})
##
## show-method for objects of class NeosJob
##
setMethod("show", "NeosJob", function(object){
  cat("\n")
  cat(paste("The job number is:", object@jobnumber, "\n"))
  cat(paste("The pass word is :", object@password, "\n"))
  cat("\n")  
})
##
## show-method for objects of class NeosOff
##
setMethod("show", "NeosOff", function(object){
  title <- paste("# The new offset is:", object@offset, "#", sep=" ")
  row <- paste(rep("#", nchar(title)), collapse="")
  cat("\n")
  cat(object@ans)
  cat("\n")
  cat(row, "\n")
  cat(title, "\n")
  cat(row, "\n")  
  cat("\n")
})
