Td <- function(d) {
  
  trt <- max(d)
  des <- as.vector(t(d))
  dtrt <- diag(trt)
  dtrt[des, ]
} 
