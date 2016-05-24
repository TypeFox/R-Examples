"lmom.ub" <-
function(x) {
  n <- length(x)
  if(n == 1) stop("use mean() for data with one value")
  if(n < 5)  stop("a minimum of 5 data values are required
                   because 5 lmoments are to be computed")
  if(length(unique(x))  == 1) stop("all values are equal--lmoments can not be computed")
  x <- sort(x)
  L1 = 0; L2 = 0; L3 = 0; L4 = 0; L5 = 0
  for(i in seq(1,n)) {
    CL1 <- i-1
    CL2 <- CL1 * (i-1-1) / 2
    CL3 <- CL2 * (i-1-2) / 3
    CL4 <- CL3 * (i-1-3) / 4
    CR1 <- n-i
    CR2 <- CR1 * (n-i-1) / 2
    CR3 <- CR2 * (n-i-2) / 3
    CR4 <- CR3 * (n-i-3) / 4     
    L1  <- L1 + x[i]
    L2  <- L2 + x[i] * (CL1 - CR1)
    L3  <- L3 + x[i] * (CL2 - 2*CL1*CR1 + CR2)
    L4  <- L4 + x[i] * (CL3 - 3*CL2*CR1 + 3*CL1*CR2 - CR3)
    L5  <- L5 + x[i] * (CL4 - 4*CL3*CR1 + 6*CL2*CR2 - 4*CL1*CR3 + CR4)    
  }
  
  C1 <- n
  C2 <- C1 * (n-1) / 2
  C3 <- C2 * (n-2) / 3
  C4 <- C3 * (n-3) / 4
  C5 <- C4 * (n-4) / 5
  L1 <- L1 / C1
  L2 <- L2 / C2 / 2
  L3 <- L3 / C3 / 3
  L4 <- L4 / C4 / 4
  L5 <- L5 / C5 / 5
  z <- list(L1 = L1, L2 = L2, TAU3 = L3/L2, TAU4 = L4/L2, TAU5 = L5/L2,
            LCV = L2/L1, L3 = L3, L4 = L4, L5=L5,
            source = "lmom.ub"
            )
  return(z)
}

