read.essentia <- function(file) {
  pipedline <- pipe(sprintf("bash %s",file), open = "r")
  colspec <- TRUE
  if (grepl("-notitle", file)) {
    colspec <- FALSE
  }
  t2 <- read.csv(pipedline, header = colspec, sep = ",", 
                 quote = "\"'", comment.char = "#", blank.lines.skip = FALSE, 
                 allowEscapes = TRUE, skip = 0)
  rm(colspec)
  close(pipedline)
  return(t2)
}