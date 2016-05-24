###################################################################
# cjgb
# 20120811
# Reads the INE deaths microdata in its 2011 (tested) version
#	It could possibly also read microdata files for other years
###################################################################

defun2011 <- function( file ){
  
  file.column  <- create.spss.column(system.file( "metadata", "defun_2011_mdat1.txt", package = "MicroDatosEs" ), 
                                     system.file( "metadata", "defun_2011_mdat2.txt", package = "MicroDatosEs" ), fileEncoding = "UTF-8")
  file.var     <- create.spss.var(system.file( "metadata", "defun_2011_mdat1.txt", package = "MicroDatosEs" ), fileEncoding = "UTF-8")
  file.vals    <- create.spss.vals(system.file( "metadata", "defun_2011_mdat2.txt", package = "MicroDatosEs" ), fileEncoding = "UTF-8")
  file.missing <- system.file( "metadata", "defun_2011_mdat3.txt", package = "MicroDatosEs" )
  
  defun2011 <- spss.fixed.file( 
                  file = file,
                  columns.file = file.column,
                  varlab.file = file.var,
                  missval.file = file.missing,
                  codes.file  = file.vals )
  
  as.data.set(defun2011)
}

