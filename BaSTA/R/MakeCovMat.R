MakeCovMat <-
function(x, data){
  if (colnames(data)[1] %in% c("id", "Id", "ID")) {
    ID <- data[, 1]
    data <- data[, -1]
  } else {
    ID <- rownames(data)
    if (is.null(ID)) {
      ID <- 1:nrow(data)
    }
  }
  if(missing(x)){
    x <- 1:ncol(data)
  }
  if (is.character(x) & length(grep("~", x)) > 0) x <- as.formula(x)
  if (is.numeric(x) | is.character(x)) {
    if (is.numeric(x)) {
      if (!all(x %in% 1:ncol(data))) {
        stop("Some arguments in 'x' do not match the column numbers in data ",
             "frame 'data'.\n", call. = FALSE)
      } else {
        xname <- colnames(data)[x]
      }
    } else if (is.character(x)) {
      if (!all(x %in% colnames(data))) {
        stop("Some arguments in 'x' do not match the column names in data ",
             "frame 'data'.\n", call. = FALSE)
      } else {
        xname <- x
      }
    } 
    classes <- sapply(xname, function(j) class(data[, j]))
    if ("factor" %in% classes) {
      facts <- which(classes == "factor")
      if (length(facts) == 1){
        xfac <- xname[classes == "factor"]
      } else {
        xfac <- paste(xname[classes == "factor"], collapse = ":")
      } 
    } else {
      xfac <- NULL
    }
    nums <- which(classes == "numeric")
    if (length(nums) == 1) {
      xnum <- xname[classes == "numeric"]
    } else if(length(nums) > 1) {
      xnum <- paste(xname[classes == "numeric"], collapse = "+")
    } else {
      xnum <- NULL
    }
    covs <- as.formula(paste("~", paste(c(xfac, xnum), 
                               collapse = "+"), "- 1"))
  } else if (class(x) == "formula") {
    xcov <- labels(terms(x))
    xcov <- unique(unlist(strsplit(xcov, split = ":")))
    if (!all(xcov %in% colnames(data))) {
      stop("Some variables in formula 'x' do not match the ",
           "column names in data.\n", call. = FALSE)
    } else {
      if (length(grep("1", x)) > 0) {
        covs <- x
      } else {
        covs <- as.formula(paste("~", paste(x, " - 1")[2], sep = ""))
      }
    }
  }
  covmat <- data.frame(ID, model.matrix(covs, data = data))
  #mode(covmat) <- "numeric"
  return(covmat)
}

