hist.influenceSSN <-
function(x, VariableName = "_resid_", xlab, ...)
{
  if(missing(xlab)) xlab = paste("VariableName = ", VariableName)

  if(class(x) != "influenceSSN") return("Not a influenceSSN object")
  hist(x$ssn.object@obspoints@SSNPoints[[1]]@point.data[,VariableName], xlab=xlab, main="", ...)

}
