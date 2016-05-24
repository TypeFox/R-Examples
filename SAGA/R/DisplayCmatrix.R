DisplayCmatrix <- function(table = "XY"){
  if(table == "esd"){
    Cmatrix <- read.csv(file = system.file("cmatrix.esd.csv", package = "SAGA"), row.names=1)[, -1]
  } else if (table == "XY"){
    Cmatrix <- read.csv(file = system.file("cmatrix.xy.csv", package = "SAGA"), row.names=1)[, -1]
  } else if (table == "ZW"){
    Cmatrix <- read.csv(file = system.file("cmatrix.zw.csv", package = "SAGA"), row.names=1)[, -1]
  } else {
    stop(cat("No genetic contribution table was found matching your argument"))
  }
  return(Cmatrix)
}
