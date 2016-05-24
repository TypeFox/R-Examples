# a temporary function to switch to the old legacy world map (1990)
world.legacy <- function(wl=FALSE){
  if(wl) Sys.setenv("R_MAP_DATA_DIR_WORLD" =
                  paste(Sys.getenv("R_MAP_DATA_DIR"),"legacy_", sep="/"))
  else Sys.setenv("R_MAP_DATA_DIR_WORLD" = Sys.getenv("R_MAP_DATA_DIR"))
  .Deprecated('map("legacy_world",...)')
}
