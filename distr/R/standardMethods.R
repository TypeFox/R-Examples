################################################################################

## standardMethods

## generates standard methods for an existing class

## args:
##
## class         -    string containing the class names
## writetofile   -    logical value: TRUE, if output is to be written on file
## directory     -    target directory for the file to be created

## value:
##
## none

## Details:
##
## If output, too, is to be written to a file, standardMethods generates a file
## <classname>_StandardMethods.txt in the given directory

################################################################################

standardMethods <- function(class, writetofile = FALSE, directory){
  ClassRep <- getClassDef(class)
  SlotNames <- names(getSlots(ClassRep))
  nrSlots <- length(SlotNames)
  
  if(!nrSlots) return()
  
  part0 <- "if(!isGeneric(\""
  part1 <- "\")) setGeneric(\""
  part2 <- "\", function(object) standardGeneric(\""
  part3 <- "\"))"
  
  for(i in 1:nrSlots){
    string <- paste(part0, SlotNames[i], part1, SlotNames[i],
                    part2, SlotNames[i], part3, "\n", sep = "")
    if(writetofile) cat(string, 
                        file = paste(directory, class, 
                                     "_StandardMethods.txt", sep=""), 
                        append = FALSE)
    cat(string, sep = "")
  }    
  
  part1 <- "setMethod(\""
  part2 <- "\", \""
  part3 <- "\", function(object) object@"
  part4 <- ")"
  
  for(i in 1:nrSlots){
    string <- paste(part1, SlotNames[i], part2, class, 
                    part3, SlotNames[i], part4, "\n", sep = "")
    if(writetofile) cat(string, 
                        file = paste(directory, class, 
                                     "_StandardMethods.txt", sep=""), 
                        append = TRUE)
    cat(string, sep = "")
  }    
  
  part0 <- "if(!isGeneric(\""
  part1 <- "<-\")) setGeneric(\""
  part2 <- "<-\", function(object, value) standardGeneric(\""
  part3 <- "<-\"))"
  
  for(i in 1:nrSlots){
    string <- paste(part0, SlotNames[i], part1, SlotNames[i], 
                    part2, SlotNames[i], part3, "\n", sep = "")
    if(writetofile) cat(string, 
                        file = paste(directory, class, 
                                     "_StandardMethods.txt", sep=""), 
                        append = TRUE)
    cat(string, sep = "")
  }    
  
  part1 <- "setReplaceMethod(\""
  part2 <- "\", \""
  part3 <- "\", function(object, value){ object@"
  part4 <- " <- value; object})"
  
  for(i in 1:nrSlots){
    string <- paste(part1, SlotNames[i], part2, class, 
                    part3, SlotNames[i], part4, "\n", sep = "")
    if(writetofile) cat(string, 
                        file = paste(directory, class, 
                                     "_StandardMethods.txt", sep=""), 
                        append = TRUE)
    cat(string, sep = "")
  }    
}            



## Example

##setClass("testclass", representation(a = "numeric", b = "character"))
##standardMethods("testclass")
##directory = "C:/TMP/"  ## directory C:\TMP must exist ...
##standardMethods("testclass", writetofile = TRUE, directory = directory)
