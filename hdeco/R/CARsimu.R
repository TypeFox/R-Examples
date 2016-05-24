"CARsimu" <-
function(level=5, row1=0.2499999, row2=0, col1=0.2499999, col2=0, rc1=0, cr1=0) {

  #-----------------------------------------------------------------------------------
  #
  #  FUNCTION:      CARsimu()
  #  DATE:          17 December, 2004
  #  AUTHOR:        SANDOR KABOS, FERENC (FERKO) CSILLAG, TARMO REMMEL
  #  CALLS:         none
  #  DESCRIPTION:   Simulates a 2^level-by-2^level real-valued landscape 
  #                 based on  the FFT algorithm and the spectral (or autocorrelation) 
  #                 theorem using first and second neighbours in N-S, E-W, 
  #                 NW-SE and NE-SW directions.
  #
  #  ARGUMENTS:     level: controls the size of the data set (2^level-by-2^level) 
  #                 row1: first-neighbour E-W autocorrelation parameter
  #                 col1: first-neighbour N-S autocorrelation parameter
  #                 row2: second-neighbour E-W autocorrelation parameter
  #                 col2: second-neighbour N-S autocorrelation parameter
  #                 rc1: first-neighbour NW-SE autocorrelation parameter
  #                 cr1: first-neighbour NE-SW autocorrelation parameter
  #                 rajz: if TRUE, an image of the otuput array is drawn on 
  #                       the graphic device
  #
  #  DETAILS:       The sum of all six autocorrelation parameters MUST be
  #                 less than 0.5. Their effect is highly non-linear, therefore
  #                 there is a marked difference between 0.499 and 0.499999.
  #                 This implementation corresponds to the 
  #                 homogeneous (stationary) conditional
  #                 autoregressive (CAR) model which has theoretical
  #                 and practical advantages over other schemes.
  #                 First-neighbour (rook's case) effects, in our experience, 
  #                 can model a wide variety of cases very well (and then 
  #                 row2=0, col2=0, cr1=0, rc=0.
  #                 Usual values for level are between 4 and 10.
  #
  #  VALUE:         Returns a 2^level-by-2^level matrix; if rajz is TRUE, also
  #                 produces an image on the graphic display
  #
  #  HISTROY:       This function was originally called ujabki()
  #
  #  REFERENCES:    Details and further technical references regarding this 
  #                 simulator can be found in section 2.2 of Remmel and Csillag (2003)
  #                 Remmel, T.K. and F. Csillag. 2003. When are two landscape 
  #                 pattern indices significantly different? 
  #                 Journal of Geographical Systems 5(4):331-351.
  #
  #  USAGE:         outarray <- CARsimu()
  #-----------------------------------------------------------------------------------


  # STEP-1: SETTING UP THE TOEPLITZ MATRIX

  M <- 2^level
  NN <- M^2
  FNN <- NN/2
  I1 <- array(rep(0, NN))
  I1[1] <- 1
  NR1 <- array(rep(0, NN))
  NR1[2] <- 1
  NR1[NN] <- 1
  NC1 <- array(rep(0, NN))
  NC1[M + 1] <- 1
  NC1[NN - M + 1] <- 1
  NR2 <- array(rep(0, NN))
  NR2[3] <- 1
  NR2[NN - 1] <- 1
  NC2 <- array(rep(0, NN))
  NC2[2 * M + 1] <- 1
  NC2[NN - 2 * M + 1] <- 1
  NRC1 <- array(rep(0, NN))
  NRC1[M + 2] <- 1
  NRC1[NN - M] <- 1
  NCR1 <- array(rep(0, NN))
  NCR1[M] <- 1
  NCR1[NN - M + 2] <- 1
  k <- I1 - row1 * NR1 - col1 * NC1 - row2 * NR2 - col2 * NC2 - rc1 * NRC1 - cr1 * NCR1
  
  # USING THE FFT
  d2 <- fft(k)
  qse <- sqrt(Re(d2))
  xe <- rnorm(NN)
  x <- fft(xe)
  for(i in 1:NN) {
    x[i] <- x[i]/qse[i]
  }

  y <- (fft(x, inv = T))/M
  yre <- Re(y)
  A <- matrix(yre, byrow = T, nrow = M, ncol = M)

  A <- A/sum(A*A)

  # ADD TITLE ATTRIBUTE 'cim' TO IDENTIFY PARAMETERS
  attr(A, "cim") <- paste("R1:", row1, " R2:", row2, " C1:", col1, " C2:", col2, " RC1:", rc1, " CR1:",cr1, sep="")  

  return(A)

}

