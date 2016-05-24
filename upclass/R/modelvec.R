.packageName <- 'upclass'

modelvec <- function (d = 1) 
{
  if (d == 1) {
    c("E", "V")
  }
  else {
    c("EII", "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "EEV", 
      "VEV", "VVV")
  }
}
