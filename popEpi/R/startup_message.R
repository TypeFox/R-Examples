
.onAttach <- function(...) {
  if (interactive()) {
    packageStartupMessage("Using popEpi. See news(package='popEpi') for changes.")
    packageStartupMessage("popEpi's appropriate data outputs are in data.table (enhanced data.frame) format by default;")
    packageStartupMessage("see ?popEpi for changing this.")
    
  }
  options("popEpi.datatable" = TRUE)
  
}




