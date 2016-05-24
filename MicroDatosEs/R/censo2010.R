###################################################################
# cjgb
# 20140428
# Reads the microdata for the 2010 Spanish Census
###################################################################

censo2010 <- function(census.file, columns = NULL, summary = TRUE){
  
  file.column  <- create.spss.column(system.file( "metadata", "censo_2010_mdat1.txt", package = "MicroDatosEs" ), 
                                     system.file( "metadata", "censo_2010_mdat2.txt", package = "MicroDatosEs" ), fileEncoding = "UTF-8")
  file.var     <- create.spss.var(system.file( "metadata", "censo_2010_mdat1.txt", package = "MicroDatosEs" ), fileEncoding = "UTF-8")
  file.vals    <- create.spss.vals(system.file( "metadata", "censo_2010_mdat2.txt", package = "MicroDatosEs" ), fileEncoding = "UTF-8")
  file.missing <- system.file( "metadata", "censo_2010_mdat3.txt", package = "MicroDatosEs" )
  
  x <- spss.fixed.file( 
              file = census.file,
              columns.file = file.column,
              varlab.file = file.var,
              missval.file = file.missing,
              codes.file  = file.vals )
  
  if(summary){
    print(summary(x))
    invisible(NULL)
  }
  
  if(! is.null(columns)){
    if(! is.character(columns))
      stop("Parameter columns needs to be a character string (column names)")
    
    columns <- names(x) %in% columns
    return(subset(x, subset = columns))
  }

  as.data.set(x)
}

