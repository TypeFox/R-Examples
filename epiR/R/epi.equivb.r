epi.equivb <- function(treat, control, delta, n, r = 1, power, alpha){

  # alpha <- (1 - conf.level)
  beta <- (1 - power)

  z.alpha <- qnorm(1 - alpha, mean = 0, sd = 1)
  z.beta <- qnorm(1 - beta / 2, mean = 0, sd = 1)

  if (!is.na(treat) & !is.na(control) & !is.na(delta) & !is.na(power) & is.na(n)) {
     # http://powerandsamplesize.com/Calculators/Test-1-Proportion/1-Sample-Equivalence:
     n.control <- (((treat * (1 - treat)) / r) + (control * (1 - control))) * ((z.alpha + z.beta) / (abs(treat - control) - delta))^2
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
     
     z <- (abs(treat - control) - delta) / sqrt((treat * (1 - treat) / n.treat) + (control * (1 - control) / n.control))
     power = 2 * (pnorm(z - qnorm(1 - alpha)) + pnorm(-z - qnorm(1 - alpha))) - 1
     
     rval <- list(n.treat = n.treat, n.control = n.control, n.total = n.treat + n.control, power = power)
  }
  rval
}  

# Chow S, Shao J, Wang H. 2008. Sample Size Calculations in Clinical Research. 2nd Ed. Chapman & Hall/CRC Biostatistics Series. page 89

# epi.equivb(treat = 0.65, control = 0.85, delta = 0.05, n = NA, power = 0.80, r = 1, alpha = 0.05)
# n.treat = 136, n.control = 136, n.total = 272
# Agrees with http://powerandsamplesize.com/Calculators/Compare-2-Proportions/2-Sample-Equivalence

# epi.equivb(treat = 0.65, control = 0.85, delta = 0.05, n = NA, power = 0.80, r = 1, alpha = 0.05)
# n.treat = 136, n.control = 136, n.total = 272
# Agrees with https://www.sealedenvelope.com/power/binary-equivalence/

# epi.equivb(treat = 0.65, control = 0.85, delta = 0.05, n = 200, power = NA, r = 1, alpha = 0.05)