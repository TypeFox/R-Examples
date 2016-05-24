qqnorm.influenceSSN <-
function(y, VariableName = "_resid_", ...)
{
  qqnorm(y$ssn.object@obspoints@SSNPoints[[1]]@point.data[,VariableName], ...)

}
