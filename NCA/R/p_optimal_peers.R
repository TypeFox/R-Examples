p_optimal_peers <-
function (input, x, y) {
  ToMatrix <- function (value) {
    tmp <- value
    if (is.vector(tmp) == TRUE) {
      tmp <- t(as.matrix(tmp))
    }
    return(tmp)
  }
  
  # Find  the optimal efficiencies for each data point
  K   <- Benchmarking::eff(input)
  K1  <- cbind (K, x, y) # add X and Y coordinates of each data point to the matrix
  NK1 <- length(K1[,2]) # find the number of rows
  K2  <- cbind (K1, c(1:NK1)) #add a number for each data point
  K3  <- K2[, c(4, 1, 2, 3)] # place these numbers as first colomn
  K4  <- subset(K3, K3[,2]> 0.999999 & K3[,2]<1.0000001) #select the peers
  
  # Select peer with highest Y  for  peers  with same X
  # sort the peer matrix according to decreasing Y- value .
  K4a <- ToMatrix(K4[sort.list(-K4[,4]), ])
  
  # sort the peer matrix according to increasing X- value.
  # Now for each set of rows with the same X,
  # the row with the highest Y values will be on the top of this set of rows.
  K4b <- ToMatrix(K4a[sort.list(K4a[,3]), ])
  
  # select for each set the row with the highest Y
  # delete duplications for the same X
  # the first row of set of rows with same X is maintained,
  # and rows with unique X are maintained
  K4c <- ToMatrix(K4b[!duplicated(K4b[,3]),])
  
  # sort the corrected peer matrix according to increasing X- value
  K4d <- ToMatrix(K4c[sort.list(K4c[,3]), ])
  
  # Select peer with lowest X  for  peers  with same Y
  # sort the peer matrix according to increasing X- value .
  K4e <- ToMatrix(K4d[sort.list(K4d[,3]), ])
  
  # sort the peer matrix according to decreasing Y- value.
  # Now for each set of rows with the same Y,
  # the row with the lowest X values will be on the top of this set of rows.
  K4f <- ToMatrix(K4e[sort.list(K4e[,3]), ])
  
  # select for each set the row with the lowest X
  # delete duplications for the same Y
  # the first row of set of rows with same Y is maintained,
  # and rows with unique Y are maintained
  K4g <- ToMatrix(K4f[!duplicated(K4f[,4]),])

  peers <- K4g[sort.list(K4g[,3]), ]

  return ( peers[, c(3, 4)] )
}