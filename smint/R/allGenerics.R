##****************************************************************************
## Despite of its name, this file DOES NOT contain ALL generics.
## 
## Due to some problems, the 'sampleIn' generic is defined in the file
## Grid.R.
##
##****************************************************************************

##=======================================================================
## Extract the official name of (USUALLY ONE-DIMENSIONAL) kernel. 
##=======================================================================
## if (!isGeneric("dim")) {
##   setGeneric("dim", function(x) standardGeneric("dim"))
## }

## if (!isGeneric("dimnames")) {
##    setGeneric("dimnames", function(x) standardGeneric("dimnames"))
## }

## if (!isGeneric("levels")) {
##    setGeneric("levels", function(x) standardGeneric("levels"))
## }

## if (!isGeneric("nlevels")) {
##    setGeneric("nlevels", function(x) standardGeneric("nlevels"))
## }

## if (!isGeneric("index")) {
##    setGeneric("index", function(object) standardGeneric("index"))
## }

if (!isGeneric("as.data.frame")) {
   setGeneric("as.data.frame",
              function(x, row.names = NULL, optional = FALSE, ...) {
                 standardGeneric("as.data.frame")
              })
}
if (!isGeneric("as.matrix")) {
   setGeneric("as.matrix",
              function(x, ...) {
                 standardGeneric("as.matrix")
              })
}
##****************************************************************************
if (!isGeneric("aperm")) {
  setGeneric("aperm",
             function(a, perm, ...){
               standardGeneric("aperm")
             })
}

