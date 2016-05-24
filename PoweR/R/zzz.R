.onLoad <- function(lib, pkg){

  library.dynam("PoweR", pkg, lib)

}

.onAttach <- function(lib, pkg){

  env <- as.environment("package:PoweR")
  
  source(system.file(package="PoweR","laws","densities.R"),env)
  source(system.file(package="PoweR","laws","moments.R"),env)
#  for (file in list.files(system.file(package="PoweR","printTemplates"),pattern="^.*\\.R$")) source(system.file(package="PoweR","printTemplates",file),env)
  
}
