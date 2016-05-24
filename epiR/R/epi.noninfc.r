epi.noninfc <- function(treat, control, sd, delta, n, r = 1, power, alpha){
  
  # alpha <- (1 - conf.level)
  beta <- (1 - power)

  z.alpha <- qnorm(1 - alpha, mean = 0, sd = 1)
  z.beta <- qnorm(1 - beta, mean = 0, sd = 1)

  if (!is.na(treat) & !is.na(control) & !is.na(delta) & !is.na(power) & is.na(n)) {
      # http://powerandsamplesize.com/Calculators/Compare-2-Means/2-Sample-Non-Inferiority-or-Superiority:
      n.control <- (1 + 1 / r) * (sd * (z.alpha + z.beta) / (treat - control - delta))^2
      n.treat <- n.control * r
     
      # r = n.treat / n.control
      n.control <- ceiling(n.control)
      n.treat <- ceiling(n.treat)
      
      rval <- list(n.treat = n.treat, n.control = n.control, n.total = n.treat + n.control)
  }
  
  if (!is.na(treat) & !is.na(control) & !is.na(delta) & !is.na(n) & is.na(power) & !is.na(r) & !is.na(alpha)) {
      # Work out the number of subjects in the control group. r equals the number in the treatment group divided by 
      # the number in the control group.
      n.control <- ceiling(1 / (r + 1) * (n))
      n.treat <- n - n.control

      
      z <- (treat - control - delta) / (sd * sqrt((1 + 1 / r) / n.control))
      power <- pnorm(z - z.alpha) + pnorm(-z - z.alpha)
      
      rval <- list(n.treat = n.treat, n.control = n.control, n.total = n.treat + n.control, power = power)
  }
  rval
}  

# Julious SA. Sample sizes for clinical trials with Normal data. Statist. Med. 2004; 23:1921-1986.
# Chow S, Shao J, Wang H. 2008. Sample Size Calculations in Clinical Research. 2nd Ed. Chapman & Hall/CRC Biostatistics Series. page 61.

# epi.noninfc(treat = 5, control = 5, sd = 10, delta = 5, n = NA, power = 0.80, r = 1, alpha = 0.05)
# n.treat = 50, n.control = 50, n.total = 100
# Agrees with http://powerandsamplesize.com/Calculators/Compare-2-Means/2-Sample-Non-Inferiority-or-Superiority
