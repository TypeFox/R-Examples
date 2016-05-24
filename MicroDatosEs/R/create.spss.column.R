###################################################################
# cjgb
# 20120305
# Transforms metadata file into data list spss file
###################################################################

create.spss.column <- function(file.start.end, file.type, ...){
  mdat.1 <- read.table(file.start.end, header = T, sep = "\t",...)

  # detecting variables whose codes are not numeric
  tmp <- read.table(file.type, header = T, sep = "\t",...)
  tmp <- tmp[tmp$tipo == "D", c("var", "llave")]
  tmp <- tmp[!grepl("^[0-9]*$", tmp$llave),]
  tmp <- as.character(unique(tmp$var))

  mdat.1$is.alpha <- ""
  mdat.1$is.alpha[mdat.1$var %in% tmp] <- "(a)"
  
  spss.column.file <- tempfile(pattern = "file", tmpdir = tempdir(), fileext = ".sps")

  cat('DATA LIST /\n ',
      paste(mdat.1$var, ' ', mdat.1$start, '-', mdat.1$end, " ", mdat.1$is.alpha, '\n ', sep = ''),
      '. \n',
      file = spss.column.file)
  
  spss.column.file
}
