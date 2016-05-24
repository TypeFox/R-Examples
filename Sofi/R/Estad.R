Estadistica <- function(Nivel){
  if (Nivel==1) {shinyAppDir(system.file("Estadist/Cal_Dis", package="Sofi"))}
  else if (Nivel==2) {shinyAppDir(system.file("Estadist/Distrib", package="Sofi"))}
  else if (Nivel==3) {shinyAppDir(system.file("Estadist/Estas", package="Sofi"))}
  else if (Nivel==4) {shinyAppDir(system.file("Estadist/ggMarginal", package="Sofi"))}
}
