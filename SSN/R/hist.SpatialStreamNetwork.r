hist.SpatialStreamNetwork <-
function(x, VariableName, main, xlab, ...)
{

  if(class(x) != "SpatialStreamNetwork") stop("Not a SpatialStreamNetwork object")

  if(missing(VariableName)) {
      VariableNames <- colnames(x@obspoints@SSNPoints[[1]]@point.data)
      Classes <- sapply(x@obspoints@SSNPoints[[1]]@point.data,class)
      ind <- which(Classes=="numeric")
      VariableNames <- VariableNames[ind]
      Classes <- Classes[ind]
      drop <- grep("albers",VariableNames)
      if(length(drop)) VariableNames <- VariableNames[-drop]
      if(length(VariableNames)==0) stop("Choose a numeric variable name")
      VariableName <- VariableNames[1]
      msg <- paste("Please specify VariableName. First numeric variable (",VariableName,") chosen by default.\n",sep="")
      cat(msg)
  }

  if(missing(main)) main <- ""
  if(missing(xlab)) xlab <- paste("VariableName = ",VariableName)

  hist(x@obspoints@SSNPoints[[1]]@point.data[,VariableName],
       xlab = xlab, main=main, ...)

}
