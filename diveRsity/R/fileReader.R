################################################################################
# Master file reader
################################################################################
fileReader <- function(infile){
  if (typeof(infile) == "list") {
    return(infile)
  } else if (typeof(infile) == "character") {
    flForm <- strsplit(infile, split = "\\.")[[1]]
    ext <- flForm[[length(flForm)]]
    if (ext == "arp") {
      convRes <- arp2gen(infile)
      if (!is.null(convRes)) {
        cat("Arlequin file converted to genepop format! \n")
        infile <- paste(flForm[1], ".gen", sep = "")
      } else {
        infile <- paste(flForm[1], ".gen", sep = "")
      }
    }
    dat <- scan(infile, sep = "\n", what = "character", quiet = TRUE)
    if(length(strsplit(dat[4], split = "\\s+")[[1]][-1]) > 1){
      locs <- strsplit(dat[2], split = "\\s+")[[1]]
      if(length(locs != 1)){
        locs <- strsplit(dat[2], split = ",")[[1]]
      }
      locs <- as.character(sapply(locs, function(x){
        x <- strsplit(x, split = "")[[1]]
        if(is.element(",", x)){
          x <- x[-(which(x == ","))]
        }
        return(paste(x, collapse = ""))
      }))
      dat <- c(dat[1], locs, dat[-(1:2)])
    }
    
    
    popLoc <- grep("^([[:space:]]*)pop([[:space:]]*)$", tolower(dat))
    no_col <- popLoc[1] - 1
    if (popLoc[1] == 3) {
      locs <- unlist(strsplit(dat[2], split = c("\\,", "\\s+")))
      dat <- c(dat[1], locs, dat[3:length(dat)])
    }
    popLoc <- grep("^([[:space:]]*)pop([[:space:]]*)$", tolower(dat))
    no_col <- popLoc[1] - 1
    dat1 <- sapply(dat, function(x) {
      x <- unlist(strsplit(x, split = "\\s+"))
      if (is.element("", x)) {
        x <- x[-(which(x == ""))]
      }
      if (is.element(",", x)) {
        x <- x[-(which(x == ","))]
      }
      if (length(x) != 1 && length(x) != no_col) {
        x <- paste(x, collapse = "")
      }
      if (length(x) < no_col) {
        tabs <- paste(rep(NA, (no_col - length(x))), 
                      sep = "\t", collapse = "\t")
        line <- paste(x, tabs, sep = "\t")
        line <- unlist(strsplit(line, split = "\t"))
        return(line)
      } else {
        return(x)
      }
    })
  }
  out <- as.data.frame(t(dat1))
  rownames(out) <- NULL
  return(out)
}
################################################################################
# END
################################################################################