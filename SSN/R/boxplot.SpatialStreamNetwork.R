boxplot.SpatialStreamNetwork <-
function(x, variable, ...)
{
    if(class(x) != "SpatialStreamNetwork") stop("Not a SpatialStreamNetwork object")

    if(missing(variable)) {
      VariableNames <- colnames(x@obspoints@SSNPoints[[1]]@point.data)
      Classes <- sapply(x@obspoints@SSNPoints[[1]]@point.data,class)
      ind <- which(Classes=="numeric")
      VariableNames <- VariableNames[ind]
      Classes <- Classes[ind]
      drop <- grep("albers",VariableNames)
      VariableNames <- VariableNames[-drop]
      if(length(VariableNames)==0) stop("Choose a numeric variable name")
      variable <- VariableNames[1]
      msg <- paste("Please specify variable. First numeric variable (",variable,") chosen by default.\n",sep="")
      cat(msg)
  }

    data1 <- x@obspoints@SSNPoints[[1]]@point.data
    if(class(variable) == "formula") return(boxplot(variable, data = data1, ...))
    if(is.character(variable)) return(boxplot(data1[,variable], ...))
    return("Error: First argument must be character name or formula")
}
