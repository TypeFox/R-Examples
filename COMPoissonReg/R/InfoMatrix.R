InfoMatrix <- function(I.beta, I.betanu, I.nu){

# Put the components/blocks that comprise the full information matrix together
  I1 <- cbind(I.beta,I.betanu)
  I2 <- cbind(t(I.betanu),I.nu)

# Create the information matrix
  I <- rbind(I1,I2)

return(I)
}

