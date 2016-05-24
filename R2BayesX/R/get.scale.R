get.scale <- function(files, dir)
{
  files <- files[!grepl("Ii11iI", files, fixed = TRUE)]
  var <- NULL
  if(any(grep("scale.res", files))) {
    sample <- NULL
    if(length(var <- grep("scale", files, value = TRUE))) {
      sc <- grep("sample", var)
      if(any(sc)) {
        sample <- grep("sample", var, value = TRUE)
        sample <- read.table(paste(dir, "/", sample, sep = ""), header = TRUE)
        sample$intnr <- NULL
        sample <- as.numeric(as.matrix(sample))
        id <- 1L:length(var)
        var <- var[id != sc]
      } else var <- var[1L]
    }
    var <- paste(dir, "/", var, sep = "")
    var <- read.table(var, header = TRUE)
    sn <- colnames(var)
    var <- as.matrix(var)
    if(length(sn) < 2L) {
      if(sn=="scale")
        sn <- "Scale"
      else
        sn <- "Sigma2"
      rownames(var) <- colnames(var) <- sn
    } else rownames(var) <- "Sigma2"
    attr(var,"sample") <- sample
  }

  return(var)
}

