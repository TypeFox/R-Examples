xport.obs.header <- function()
{
  .C("fill_obs_header", PACKAGE="SASxport")
  .Call("getRawBuffer", PACKAGE="SASxport")
}

