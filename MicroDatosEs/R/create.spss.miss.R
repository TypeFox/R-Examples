###################################################################
# cjgb
# 20120305
# Transforms metadata file into data list spss file
###################################################################

create.spss.miss <- function(file,...){
  return( file )

  # completar el resto...
  mdat.1 <- read.table(file, header = T, sep = "\t",...)
  
  spss.column.file <- tempfile(pattern = "file", tmpdir = tempdir(), fileext = ".sps")
  
  cat('DATA LIST /\n ',
      paste(mdat.1$var, ' ', mdat.1$start, '-', mdat.1$end,'\n ', sep = ''),
      '. \n',
      file = spss.column.file)
  
  spss.column.file
}
