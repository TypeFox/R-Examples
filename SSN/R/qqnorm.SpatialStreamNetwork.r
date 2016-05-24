qqnorm.SpatialStreamNetwork <-
function(y, VariableName = "_resid_", ...)
{
  qqnorm(y@obspoints@SSNPoints[[1]]@point.data[,VariableName], ...)

}
