# Calculate pdf
get.pdf <- function(inputData, od, J, M, param){
  tau <- get.tau(inputData)
  pdf <- rep(0, times = J * tau)
  for (j in 1:J){
    # od = "Bernoulli"
    if (od == "bern"){
      for (t in 1:tau){
        if (inputData[t] == 1){
          pdf[t + (j - 1) * tau] <- param$b[j]
          }
        else {
          pdf[t + (j - 1) * tau] <- 1 - param$b[j]
          }
        }
      }
    # od = "Gaussian"
    if (od == "norm"){
      pdf[(1 + (j - 1) * tau):(tau + (j - 1) * tau)] <- dnorm(inputData[1:tau], mean = param$mean[j], sd = sqrt(param$var[j]))
      }
    # od = "Poisson"
    if (od == "pois"){
      pdf[(1 + (j - 1) * tau):(tau + (j - 1) * tau)] <- dpois(inputData[1:tau], lambda = param$lambda[j])
      }
    # od = "Student.t"
    if (od == "t"){
      pdf[(1 + (j - 1) * tau):(tau + (j - 1) * tau)] <- dtmod(inputData[1:tau], mu = param$mean[j], sigma = sqrt(param$var[j]), nu = param$df[j])
      }
     # od = "multivar.Gaussian"
    if (od == "mvnorm"){
      pdf[(1 + (j - 1) * tau):(tau + (j - 1) * tau)] <- dmvnorm(aperm(inputData[,1:tau]), mean = param$mean[,j], sigma = param$sigma[,,j])
      }
   }

  lower_bound <- 1e-300
  pdf[pdf < lower_bound] <- lower_bound
  return(pdf)
}
