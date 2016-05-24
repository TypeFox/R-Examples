williams <- function(trt) {
  
  if (is.numeric(trt) == FALSE || any(length(trt) != 1, trt%%1 != 0, trt < 2) == 
    TRUE) {
    stop("Number of treatments is not an integer >1.")
  }
  # Checks for appropriate values of t (must be an integer larger than 1)
  
  if (trt == 2) {
    design <- matrix(c(1, 2, 2, 1), 2, 2)
  }
  # for t=2, the williams design is matrix(c(1,2,2,1),2,2).
  if (trt > 2) {
    
    gen.row <- c(0, rep(1, trt - 1))
    a <- ifelse((3:trt)%%2 == 0, -1, 1)
    
    for (i in 3:trt) {
      gen.row[i] <- gen.row[i - 1] + a[i - 2] * (trt + 1 - i)
    }
    
    design <- rbind(gen.row, matrix(rep(0, trt * (trt - 1)), ncol = trt))
    row.names(design) <- NULL
    
    for (i in 2:trt) {
      design[i, ] <- (design[i - 1, ] + 1)%%trt
    }
    
    design <- design + 1
    
    if (trt%%2 == 1) {
      design <- rbind(design, t(apply(design, 1, rev)))
    }
    row.names(design) <- NULL
    
  }
  
  design
  
} 
