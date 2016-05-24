# --------------------------------------------
# Function for the simulated test decision
# --------------------------------------------

# Returns the exact/simulated test decision of Student's t-test
#
# Calculates for every randomization sequence the exact/simulated p.value.
#
# @param randSeq object of the class randSeq.
# @param bias object of the class bias.
# @param endp object of the class endpoint.
# 
# @return
# vector of the simulated/exact p.value of a randomization sequence.
testDec <- function(randSeq, bias, endp) {
  stopifnot(is(randSeq, "randSeq"), randSeq@K == 2,
            #is(bias, "chronBias") || is(bias, "selBias") || is(bias, "power"), 
            is(endp, "normEndp"))
  if (bias@method == "sim") {
    # calculates the bias matrix
    biasM <- getExpectation(randSeq, bias, endp)
    # matrix of the standard deviations
    sdM <- matrix(numeric(0), ncol = dim(randSeq@M)[2], nrow = dim(randSeq@M)[1])
    sdM[randSeq@M == 0] <- endp@sigma[1]
    sdM[randSeq@M == 1] <- endp@sigma[2]
    
    sapply(1:dim(randSeq@M)[1], function(i) {      
      randVar <- rnorm(length(biasM[i, ]) , mean = biasM[i, ], sd = sdM[i, ] )
      if (sum(randSeq@M[i,]) == 0 || sum(randSeq@M[i,]) == length(biasM[i, ]) ) {
        return(FALSE)
      } else {
        t.test(randVar ~ randSeq@M[i, ], var.equal = TRUE)$p.value <= bias@alpha
      }
    } 
    )
  } else if (bias@method == "exact") {
    doublyTValues(randSeq, bias, endp)
  }       
}
