setClass("dataFrame1",
  representation(row.names = "character", names = "character"),
  contains = "list")

setClass("dataFrame2", 
  representation(variables = "list",
    colNames = "character", rowNames = "character"))

setClassUnion("dataFrame",
        c("data.frame", "dataFrame1", "dataFrame2"))

## requireMethods("dimnames", "dataFrame")
## requireMethods("dim", "dataFrame")

setMethod("dim", "dataFrame1", function(x) {
    nr = length(x@row.names)
    if(nr == 0 && length(x) > 0)
      nr = length(x[[1]])
    c(nr, length(x))
})

setMethod("[", "dataFrame1", base::`[.data.frame` )#]

setValidity("dataFrame", function(object) {
    msg <- character()
    dims <- dim(object)
    if(length(dims) != 2)
      return(paste("dim() not of length 2 (got", length(dims), ")"))
    nRow <- dims[[1]]
    nVar <- dims[[2]]
    if(nVar > 0) {
        varLens <- numeric(nVar)
        for(j in seq(length = nVar))
          varLens[[j]] <- length(object[, j])
        nR2 = unique(varLens)
        if(length(nR2) > 1)
          msg <-  c(msg,
                    paste("Variables (columns) not of equal length: ", paste(nR2, collapse = ", ")))
        if(any(is.na(match(nR2, nRow))))
          msg <- c(msg,
                   paste("Length of some columns != ", nRow, " (dim()[[1]]): ",
                         paste(nR2[is.na(match(nR2, nRow))], collapse = ", ")))
    }
    if(length(msg)> 0)
      msg
    else
      TRUE
})

try(
dd = new("dataFrame1", list(a=1:10, b=rnorm(10), c=rnorm(9)), row.names = letters[1:10])
    )

