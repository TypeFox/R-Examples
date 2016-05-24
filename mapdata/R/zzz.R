.onLoad <- function(lib, pkg) {
  Sys.setenv("R_MAPDATA_DATA_DIR"=paste(lib, pkg, "mapdata/", sep="/"))
}
