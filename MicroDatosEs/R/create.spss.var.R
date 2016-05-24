###################################################################
# cjgb
# 20120305
# Transforms metadata file into data list spss file
###################################################################

create.spss.var <- function(file,...){
  mdat.1 <- read.table(file, header = T, sep = "\t",...)
  
  spss.var.file <- tempfile(pattern = "file", tmpdir = tempdir(), fileext = ".sps")
  
  cat('VARIABLE LABELS / \n',
      paste(' /', mdat.1$var,' "', mdat.1$descr,'"\n',sep=''),
      '.',
      file = spss.var.file )
  
  spss.var.file
}