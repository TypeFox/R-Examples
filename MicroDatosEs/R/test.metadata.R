###################################################################
# cjgb
# 20120305
# Function intended for testing metadata files
###################################################################

test.metadata <- function( file, md.1, md.2, md.3, encoding = "UTF-8" ){
  
  file.column  <- create.spss.column(md.1, md.2, encoding = encoding)
  file.var     <- create.spss.var(md.1, encoding = encoding)
  file.vals    <- create.spss.vals(md.2, encoding = encoding)
  file.missing <- create.spss.miss(md.3, encoding = encoding)
  
  res <- spss.fixed.file( file = file,
    columns.file = file.column,
    varlab.file = file.var,
    missval.file = file.missing,
    codes.file  = file.vals )
  
  as.data.set(res)
}
