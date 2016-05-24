Restrictions.BinBin <- function(pi1_1_, pi1_0_, pi_1_1, pi_1_0, pi0_1_, pi_0_1) {   
  
  # no mono
  min_pi_1001 <- min(pi1_0_, pi_0_1)
  min_pi_1110 <- min(pi1_1_, pi_1_0)
  min_pi_1101 <- min(pi1_0_, pi_1_1)
  min_pi_1011 <- min(pi1_1_, pi_0_1)
  min_pi_1111 <- min(pi1_1_, pi_1_1)
  min_pi_0110 <- min(pi0_1_, pi_1_0)
  min_pi_0011 <- min(pi0_1_, pi_0_1)
  min_pi_0111 <- min(pi0_1_, pi_1_1)
  min_pi_1100 <- min(pi1_0_, pi_1_0)
  
  cat("\nAssuming no monotonicity, the following restrictions on the freely varying parameters")
  cat("\nare imposed by the data:")
  cat("\n#------------------------------------------------------------------------------------")
  cat("\npi_1001 <=", min_pi_1001)
  cat("\npi_1110 <=", min_pi_1110)
  cat("\npi_1101 <=", min_pi_1101)
  cat("\npi_1011 <=", min_pi_1011)
  cat("\npi_1111 <=", min_pi_1111)
  cat("\npi_0110 <=", min_pi_0110)
  cat("\npi_0011 <=", min_pi_0011)
  cat("\npi_0111 <=", min_pi_0111)
  cat("\npi_1100 <=", min_pi_1100)
  max_val <- min((min_pi_1001 + min_pi_1110 + min_pi_1101 + min_pi_1011 + min_pi_1111 + 
                    min_pi_0110 + min_pi_0011 + min_pi_0111 + min_pi_1100), 1)
  
  cat("\n\nSum Pi_f <= ", max_val)
  
  
  # mono T
  min_pi_1111 <- min(pi1_1_, pi_1_1)
  min_pi_0110 <- min(pi0_1_, pi_1_0)
  min_pi_0011 <- min(pi0_1_, pi_0_1)
  min_pi_0111 <- min(pi0_1_, pi_1_1)
  min_pi_1100 <- min(pi1_0_, pi_1_0)
  
  cat("\n\n\nAssuming monotonicity for T, the following restrictions on the freely varying parameters")
  cat("\nare imposed by the data:")
  cat("\n#-------------------------------------------------------------------------------------------")
  cat("\npi_1111 <=", min_pi_1111)
  cat("\npi_0110 <=", min_pi_0110)
  cat("\npi_0011 <=", min_pi_0011)
  cat("\npi_0111 <=", min_pi_0111)
  cat("\npi_1100 <=", min_pi_1100)
  max_val <- min((min_pi_1111 + min_pi_0110 + min_pi_0011 + min_pi_0111 + min_pi_1100), 1)
  
  cat("\n\nSum Pi_f <= ", max_val)
  

  # mono S
  min_pi_1001 <- min(pi1_0_, pi_0_1)
  min_pi_1101 <- min(pi1_0_, pi_1_1)
  min_pi_1111 <- min(pi1_1_, pi_1_1)
  min_pi_0111 <- min(pi0_1_, pi_1_1)
  min_pi_1100 <- min(pi1_0_, pi_1_0)
  
  cat("\n\n\nAssuming monotonicity for S, the following restrictions on the freely varying parameters")
  cat("\nare imposed by the data:")
  cat("\n#-------------------------------------------------------------------------------------------")
  cat("\npi_1001 <=", min_pi_1001)
  cat("\npi_1101 <=", min_pi_1101)
  cat("\npi_1111 <=", min_pi_1111)
  cat("\npi_0111 <=", min_pi_0111)
  cat("\npi_1100 <=", min_pi_1100)
  max_val <- min((min_pi_1001 + min_pi_1101 + min_pi_1111 + min_pi_0111 + min_pi_1100), 1)
  
  cat("\n\nSum Pi_f <= ", max_val)
  
  
  # mono S and T
  min_pi_1001 <- min(pi1_0_, pi_0_1)
  min_pi_1110 <- min(pi1_1_, pi_1_0)
  min_pi_1101 <- min(pi1_0_, pi_1_1)
  min_pi_1011 <- min(pi1_1_, pi_0_1)
  min_pi_1111 <- min(pi1_1_, pi_1_1)
  min_pi_0110 <- min(pi0_1_, pi_1_0)
  min_pi_0011 <- min(pi0_1_, pi_0_1)
  min_pi_0111 <- min(pi0_1_, pi_1_1)
  min_pi_1100 <- min(pi1_0_, pi_1_0)
  
  cat("\n\n\nAssuming monotonicity for S and T, the following restrictions on the freely varying parameters")
  cat("\nare imposed by the data:")
  cat("\n#---------------------------------------------------------------------------------------------")
  min_pi_0111 <- min(pi0_1_, pi_1_1)
  min_pi_1100 <- min(pi1_0_, pi_1_0)

  cat("\npi_0111 <=", min_pi_0111)
  cat("\npi_1100 <=", min_pi_1100)
    
  max_val <- min((min_pi_0111 + min_pi_1100), 1)
  
  cat("\n\nSum Pi_f <= ", max_val)
  
  
}