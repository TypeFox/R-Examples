#' A function to write a genepop file, where individuals are grouped into population samples based on a factor or vector.
#'
#' Kevin Keenan, 2014
#'
regroup <- function(infile = NULL, index = NULL){
  dat <- rgp(infile)
  if(length(index) != length(unlist(dat$indnms))){
    stop("The number of individuals in index does not match the input file!")
  }
  if(is.null(names(index))){
    names(index) <- gsub(",", "", unlist(dat$indnms))
  }
  #if(rm.comma){
  nms <- lapply(dat$indnms, function(x){
    gsub(",", "", x)
  })
  #} else {
  #  nms <- dat$indnms
  #}
  idxs <- lapply(names(index), function(ind){
    idx_1 <- which(sapply(nms, function(y){
      any(y == ind)
    }))
    idx_2 <- which(nms[[idx_1]] == ind)
    out <- c(idx_1, idx_2)
    names(out) <- c(ind, ind)
    return(out)
  })
  # split
  gp <- round(mean(nchar(dat$genos[[1]][!is.na(dat$genos[[1]])])))
  idxs <- tapply(idxs, index, "[")
  new_gt <- lapply(idxs, function(x){
    nms <- paste(sapply(x, function(y) names(y)[1]), ",", sep = " ")
    out <- t(sapply(x, function(y){
      ret <- apply(dat$genos[[y[1]]][y[2],,], 1, paste, collapse = "")
      if(gp == 3L) {
        ret[ret == "NANA"] <- "000000"
      } else {
        ret[ret == "NANA"] <- "0000"
      }
      return(ret)
    }))
    return(cbind(nms, out))
  })
  
  hdr <- c("regrouped genepop file", dat$locs)
  new_gt <- lapply(new_gt, function(x){
    c("POP", apply(x, 1, paste, collapse = "\t"))
  })
  flForm <- strsplit(infile, split = "\\.")[[1]]
  if(substr(infile, 1, 2) == "./"){
    flForm <- flForm[-1]
  } else if(substr(infile, 1, 3) == "../"){
    flForm <- flForm[-(1:2)]
  }
  new_gt <- paste(c(hdr, unlist(new_gt)), collapse = "\n")
  writeLines(new_gt, paste(".", flForm[1], "-regroup.gen", sep = ""))
}