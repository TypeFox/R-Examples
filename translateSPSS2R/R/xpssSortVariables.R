#'  Sorting variables
#'
#'  R implementation of the SPSS \code{SORT VARIABLES} function.
#'
#' @usage xpssSortVariables(x, by = NULL, order = "A")
#' @param x a (non-empty) data.frame or input data of class "xpssFrame".  
#' @param by atomic character with the name of the argument to sort by  
#' @param order atomic character which indicate the sort direction.
#'
#' @return Returns \code{x} sorted by the \code{by} argument. 
#' 
#' @details \code{x} can be sorted based on the following \code{by} arguments. 
#'  \itemize{
#'  \item \emph{NAME.} Sort by variable names. Primary sort is alphabetical. Trailing digits are sorted numerically 
#'  \item \emph{TYPE.} Sort variables by type (numeric or string). Sort string variables by width.
#'  \item \emph{LABEL.} Sort variables alphabetical by variable labels.
#'  \item \emph{COLUMNS.} Sort variables by column width.
#'  }
#'  Valid input for order are "Up" or "Down", respectively "A" or "D".
#'
#' @author Benjamin Piest
#' @seealso \code{\link{sort}}
#' @examples 
#'
#'data(fromXPSS)
#'
#' xpssSortVariables(fromXPSS, by = "NAME")
#' 
#' xpssSortVariables(fromXPSS, by = "NAME", order = "D")
#' 
#' xpssSortVariables(fromXPSS, by = "TYPE", order = "A")
#' 
#' xpssSortVariables(fromXPSS, by = "TYPE", order = "D")
#' 
#' xpssSortVariables(fromXPSS, by = "COLUMNS", order = "A")
#' 
#' xpssSortVariables(fromXPSS, by = "COLUMNS", order = "D")
#' 
#' @export 


xpssSortVariables <- function(x, by = NULL, order = "A"){
  
  if(order == "UP"){order="A"}
  if(order == "DOWN"){order="D"}
  
#Zahlem m?ssen noch nach spss Standard sortiert werden 2 vor 10
  
  if(by == "NAME" & order == "A"){
    x <- x[, sort(names(x), F)]
  }
  if(by == "NAME" & order == "D"){
    x <- x[, sort(names(x), T)]
  }

if(by == "TYPE" & order == "A"){
  for(i in 1:ncol(x)){
    if(class(x[[i]])=="integer"){
      x[[i]] <- as.numeric(x[[i]])
    }
    if(class(x[[i]])=="character"){
      x[[i]] <- as.factor(x[[i]])
    }
    x1 <- x[,sapply(x,class)=="numeric"]
    x2 <- x[,sapply(x,class)=="factor"]
  }
  x1 <- x1[names(sort(sapply(data.frame(sapply(x1,str_length)),max,na.rm=T),F))]
  x2 <- x2[names(sort(sapply(data.frame(sapply(x2,str_length)),max,na.rm=T),F))]
  x <- data.frame(x1,x2)
}
  
if(by == "TYPE" & order == "D"){
    for(i in 1:ncol(x)){
      if(class(x[[i]])=="integer"){
        x[[i]] <- as.numeric(x[[i]])
      }
      if(class(x[[i]])=="character"){
        x[[i]] <- as.factor(x[[i]])
      }
      x1 <- x[,sapply(x,class)=="numeric"]
      x2 <- x[,sapply(x,class)=="factor"]
    }
    x1 <- x1[names(sort(sapply(data.frame(sapply(x1,str_length)),max,na.rm=T),T))]
    x2 <- x2[names(sort(sapply(data.frame(sapply(x2,str_length)),max,na.rm=T),T))]
    x <- data.frame(x2,x1)
}


if(by == "COLUMNS" & order == "A"){
  x <- x[names(sort(sapply(data.frame(sapply(x,str_length)),max,na.rm=T),F))]
}

if(by == "COLUMNS" & order == "D"){
  x <- x[names(sort(sapply(data.frame(sapply(x,str_length)),max,na.rm=T),T))]
}


if(by == "LABEL" & order == "A"){
  var.label <- 1:length(x)
  for(i in 1:length(x)){
    var.label[i] <- lapply(x, attributes)[[i]][["variable.label"]]
  }
  var.label <- matrix(var.label,1,10)
  colnames(var.label) <- names(x)
  x <- x[names(sort(var.label[1,],F))]
}

if(by == "LABEL" & order == "D"){
  var.label <- 1:length(x)
  for(i in 1:length(x)){
    var.label[i] <- lapply(x, attributes)[[i]][["variable.label"]]
  }
  var.label <- matrix(var.label,1,10)
  colnames(var.label) <- names(x)
  x <- x[names(sort(var.label[1,],T))]
}

  return(x)
}
