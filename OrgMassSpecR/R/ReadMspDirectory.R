ReadMspDirectory <- function(directory, skip = 2, comment.char = "") {

  ReadMSP <- function(file) {

    msp <- read.table(file = file, sep = ";", skip = skip, fill = TRUE, 
                      stringsAsFactors = FALSE, comment.char = comment.char)
    
    print(paste("Reading file", file))

    msp <- msp[-ncol(msp)]

    # make tall
    values <- NULL
    for (i in 1:ncol(msp)) {
      values <- c(values, msp[, i])
    }
    tall.format <- as.vector(values, mode = "character")

    RemoveWhite <- function(x) {
      y <- strsplit(x, split = "[[:space:]]")[[1]]
      z <- y[y != ""]
      return(as.numeric(z))
    }
    results.list <- lapply(tall.format, RemoveWhite)

    result <- as.data.frame(do.call("rbind", results.list))
    names(result) <- c("mz", "intensity")
    result <- result[result$intensity != 0, ]
    ordered.result <- result[order(result$mz), ]
    row.names(ordered.result) <- 1:nrow(ordered.result)
    return(ordered.result)
  
  }

  spectra.files <- dir(directory)
  spectra.list <- vector("list")
  for (i in 1:length(spectra.files)) {
    tmp <- ReadMSP(file = paste(directory, "/", spectra.files[i], sep = ""))
    tmp.filename <- strsplit(spectra.files[i], split = "\\.")[[1]][1]
    tmp$filename <- tmp.filename
  spectra.list[[i]] <- tmp
  }

  spectra.table <- as.data.frame(do.call("rbind", spectra.list))

  return(spectra.table)

}


