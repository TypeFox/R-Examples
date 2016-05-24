recode2 <- function(x1, x2, codes){
##
## 1.  check length(x1) == length(x2)
##
  N1 <- length(x1)
  N2 <- length(x2)
  if(N1 != N2)
      stop('length(x1) =', N1, '!=', N2, '= length(x2)')
##
## 2.  is.logical x1, x2
##
  if(is.logical(x1)) l1 <- c(FALSE, TRUE) else l1 <- unique(x1)
  if(is.logical(x2)) l2 <- c(FALSE, TRUE) else l2 <- unique(x1)
##
## 3.  missing(codes)
##
  if(missing(codes))
      codes. <-  outer(l1, l2, paste, sep=":")
##
## 4.  is.null(dim(codes))
##
  n1 <- length(l1)
  n2 <- length(l2)
  if(is.null(dim(codes))) dim(codes) <- c(n1, n2)
##
## 5.  is.null(rownames(codes))
##
  if(is.null(rownames(codes))){
      if(nrow(codes) == n1) {
          rownames(codes) <- l1
      } else {
          if(nrow(codes) == max(x1)){
              rownames(codes) <- 1:max(x1)
          } else stop('rownames(codes) not provided and does not match x1')
      }
  }
  if(is.null(colnames(codes))){
      if(ncol(codes) == n2) {
          colnames(codes) <- l2
      } else {
          if(ncol(codes) == max(x2)){
              colnames(codes) <- 1:max(x2)
          } else stop('colnames(codes) not provided and does not match x2')
      }
  }
##
## 6.  Do
##
  x12 <- cbind(as.character(x1), as.character(x2))
  codes[x12]
}

