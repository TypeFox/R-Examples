###################################################################
# cjgb
# 20120305
# Reads the "Encuesta de Estructura Salarial" with 2010 format
###################################################################

ees2010 <- function(ees.file){
  
  file.column  <- create.spss.column(system.file( "metadata", "ees_mdat1.txt", package = "MicroDatosEs" ), 
                                     system.file( "metadata", "ees_mdat2.txt", package = "MicroDatosEs" ), fileEncoding = "UTF-8")
  file.var     <- create.spss.var(system.file( "metadata", "ees_mdat1.txt", package = "MicroDatosEs" ), fileEncoding = "UTF-8")
  file.vals    <- create.spss.vals(system.file( "metadata", "ees_mdat2.txt", package = "MicroDatosEs" ), fileEncoding = "UTF-8")
  file.missing <- system.file( "metadata", "ees_mdat3.txt", package = "MicroDatosEs" )
  
  ees2010 <- spss.fixed.file( 
                              file = ees.file,
                              columns.file = file.column,
                              varlab.file = file.var,
                              missval.file = file.missing,
                              codes.file  = file.vals )
  
  fix.char.items(as.data.set(ees2010))
}

