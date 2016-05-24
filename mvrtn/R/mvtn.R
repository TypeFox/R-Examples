mvtn <-
function(zmu,zsig,c,side = c("left","right"))
{
  side <- match.arg(side)
  if ( side == "left")
  {
    returned_data = .C("MVTN_Left", zmu = as.double(zmu), zsig = as.double(zsig), c = as.double(c), result = as.double(numeric(2)))
    return(returned_data$result)
  }
  else if ( side == "right")
  {
    returned_data = .C("MVTN_Right", zmu = as.double(zmu), zsig = as.double(zsig), c = as.double(c), result = as.double(numeric(2)))
    returned_data$result
  }
}
