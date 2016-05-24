.onAttach <- function(...) {
  if(!interactive()) return()
  packageStartupMessage("please cite mpoly if you use it; see citation(\"mpoly\")")  
}