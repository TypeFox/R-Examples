epi.equivc <- function(treat, control, sd, delta, n, r = 1, power, alpha){
  
  # alpha <- (1 - conf.level)
  beta <- (1 - power)

  z.alpha <- qnorm(1 - alpha, mean = 0, sd = 1)
  z.beta <- qnorm(1 - beta / 2, mean = 0, sd = 1)
  
  if (!is.na(treat) & !is.na(control) & !is.na(power) & is.na(n)) {
     # http://powerandsamplesize.com/Calculators/Compare-2-Means/2-Sample-Equality:
     n <- ((sd * (z.alpha + z.beta) / (abs(treat - control) - delta))^2)
    
     # r = n.treat / n.control
     n.control <- (1  + 1 / r) * (ceiling(n))
     n.treat <- n.control * r
     
     rval <- list(n.treat = n.treat, n.control = n.control, n.total = n.treat + n.control)
  }
  
  if (!is.na(treat) & !is.na(control) & !is.na(n) & is.na(power) & !is.na(r) & !is.na(alpha)) {
     # Work out the number of subjects in the control group. r equals the number in the treatment group divided by 
     # the number in the control group.
     n.control <- ceiling(1 / (r + 1) * (n))
     n.treat <- n - n.control
     
     z <- (abs(treat - control) - delta) / (sd * sqrt((1 / n.treat) + (1 / n.control))) 
     power <- 2 * (pnorm(z - qnorm(1 - alpha)) + pnorm(-z - qnorm(1 - alpha))) - 1    
     
     rval <- list(n.treat = n.treat, n.control = n.control, n.total = n.treat + n.control, power = power)
  }
  rval
}  

# Chow S, Shao J, Wang H. 2008. Sample Size Calculations in Clinical Research. 2nd Ed. Chapman & Hall/CRC Biostatistics Series. page 62


# epi.equivc(treat = 5, control = 4, sd = 10, delta = 5, n = NA, power = 0.80, r = 1, design = 1, alpha = 0.05)
# n.treat = 108, n.control = 108, n.total = 216
# Agrees with http://powerandsamplesize.com/Calculators/Compare-2-Means/2-Sample-Equivalence

# epi.equivc(treat = 5, control = 4, sd = 10, delta = 5, n = NA, power = 0.80, r = 1, alpha = 0.05)
# n.treat = 108, n.control = 108, n.total = 216
# Agrees with https://www.sealedenvelope.com/power/continuous-equivalence/