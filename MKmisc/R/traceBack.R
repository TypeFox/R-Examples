traceBack <- function(D, global = TRUE){
  stopifnot("stringDist" %in% class(D) | "stringSim" %in% class(D))
  TB <- attr(D, "TraceBackMatrix")
  SM <- attr(D, "ScoringMatrix")
  if(global){
    m <- nrow(TB)
    n <- ncol(TB)
  }else{
    wmax <- which.max(SM)
    if(wmax %% nrow(SM) == 0){
      m <- nrow(SM)
      n <- wmax %/% nrow(SM)
    }else{
      m <- wmax %% nrow(SM)
      n <- wmax %/% nrow(SM) + 1
    }
  }
  x <- NULL
  y <- NULL
  local <- FALSE
  while(m > 1 | n > 1){
    if(length(grep("stop", TB[m,n])) > 0){
      local <- TRUE
      break
    }
    if(length(grep("m", TB[m,n])) > 0){
      x <- c(rownames(TB)[m], x)
      y <- c(colnames(TB)[n], y)
      m <- m-1
      n <- n-1
      next
    }
    if(length(grep("d", TB[m,n]))){
      x <- c(rownames(TB)[m], x)
      y <- c("-", y)
      m <- m-1
      next      
    }
    if(length(grep("i", TB[m,n]))){
      y <- c(colnames(TB)[n], y)
      x <- c("-", x)
      n <- n-1
      next      
    }
  }
  res <- rbind(paste(x, collapse = ""), paste(y, collapse = ""))
  rownames(res) <- c("x'", "y'")
  if(global)
    colnames(res) <- c("pairwise global alignment")
  else
    colnames(res) <- c("pairwise local alignment")
  res
}
