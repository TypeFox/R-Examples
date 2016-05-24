"epi.studysize" <- function(treat, control, n, sigma, power, r = 1, design = 1, sided.test = 2, conf.level = 0.95, method = "means") {
   
   alpha.new <- (1 - conf.level) / sided.test
   z.alpha <- qnorm(1 - alpha.new, mean = 0, sd = 1)
 
 if(method == "means" & !is.na(treat) & !is.na(control) & is.na(n) & !is.na(sigma) & !is.na(power)){
  # Sample size. From Woodward p 398:
  z.beta <- qnorm(power, mean = 0, sd = 1) 
  delta <- abs(treat - control)
  n <- ((r + 1)^2 * (z.alpha + z.beta)^2 * sigma^2) / (delta^2 * r)
  
  # Account for the design effect:
  n <- n * design
  
  n.crude <- ceiling(n)
  n.treat <- ceiling(n / (r + 1)) * r
  n.control <- ceiling(n / (r + 1)) * 1
  n.total <- n.treat + n.control
  rval <- list(n.crude = n.crude, n.total = n.total, n.treat = n.treat, n.control = n.control)
  }

 else 
 if(method == "means" & !is.na(treat) & !is.na(control) & !is.na(n) & !is.na(sigma) & is.na(power)){
  # Study power. From Woodward p 401:
  delta <- abs(treat - control)
  
  # Account for the design effect:
  n <- n / design
  
  z.beta <- ((delta * sqrt(n * r)) / ((r + 1) * sigma)) - z.alpha
  power <- pnorm(z.beta, mean = 0, sd = 1)
  rval <- list(power = power)
  }

 else 
 if(method == "means" & is.na(treat) & is.na(control) & !is.na(n) & !is.na(sigma) & !is.na(power)){
  # Maximum detectable difference. From Woodward p 401:
  z.beta <- qnorm(power, mean = 0, sd = 1) 
  
  # Account for the design effect:
  n <- n / design
  
  delta <- ((r + 1) * (z.alpha + z.beta) * sigma) / (sqrt(n * r))
  rval <- list(delta = delta)
  }

 else
 if (method == "proportions" & !is.na(treat) & !is.na(control) & is.na(n) & !is.na(power)) {
  # Sample size.
  z.beta <- qnorm(power, mean = 0, sd = 1)
  # delta <- abs(treat - control)
  # n <- (1/delta^2) * ((z.alpha * sqrt(treat * (1 - treat))) + (z.beta * sqrt(control * (1 - control))))^2
  
  # From Woodward's spreadsheet. Changed 130814:
  lambda <- treat / control
  Pc <- control * (r * lambda + 1) / (r + 1)
  T1 <- (r + 1) / (r * (lambda - 1)^2 * control^2)
  T2 <- (r + 1) * Pc *(1 - Pc)
  T3 <- lambda * control * (1 - lambda * control) + r * control * (1 - control)
  n <- T1 * (z.alpha * sqrt(T2) + z.beta * sqrt(T3))^2
  
  # Account for the design effect:
  n <- n * design
  
  # n.total <- 2 * ceiling(0.5 * n)
  # rval <- list(n.total = n.total)
  n.crude <- ceiling(n)
  n.treat <- ceiling(n / (r + 1)) * r
  n.control <- ceiling(n / (r + 1)) * 1
  n.total <- n.treat + n.control
  rval <- list(n.crude = n.crude, n.total = n.total, n.treat = n.treat, n.control = n.control)
  }

 else
 if (method == "proportions" & !is.na(treat) & !is.na(control) & !is.na(n) & is.na(power)) {
  # Power.
  # Account for the design effect:
  n <- n / design

  # From Woodward's spreadsheet. Changed 130814:
  lambda <- control / treat
  Pc <- treat * (r * lambda + 1) / (r + 1)
  T1 <- ifelse(lambda >= 1, treat * (lambda - 1) * sqrt(n * r), treat * (1 - lambda) * sqrt(n * r))
  T2 <- z.alpha * (r + 1) * sqrt(Pc * (1 - Pc))
  T3 <- (r + 1) * (lambda * treat * (1 - lambda * treat) + r * treat * (1 - treat))
  z.beta <- (T1 - T2) / sqrt(T3)
  # z.beta <- ((delta * sqrt(n)) - (z.alpha * sqrt(treat * (1 - treat))))/(sqrt(control * (1 - control)))
  power <- pnorm(z.beta, mean = 0, sd = 1)
  rval <- list(power = power)
  }
  
 else 
 if (method == "proportions" & !is.na(treat) & !is.na(control) & !is.na(n) & !is.na(power)) {
  # Maximum detectable difference.
  z.beta <- qnorm(power, mean = 0, sd = 1)
  
  # Account for the design effect:
  n <- n / design
  
  delta <- 1/sqrt(n) * ((z.alpha * sqrt(treat * (1 - treat))) + (z.beta * sqrt(control * (1 - control))))
  rval <- list(delta = delta)
  }

# else 
# if(method == "proportions" & is.na(n)){
#   # Sample size estimate. From Fleiss (1981).
#   z.beta <- qnorm(power, mean = 0, sd = 1) 
#   delta <- abs(treat - control)
#   n <- (z.alpha + z.beta)^2 * (((treat * (1 - treat)) + (control * (1 - control))) / delta^2) + (2 / delta) + 2
#   n <- ceiling(2 * n)
#   rval <- list(n = n)
#      }                                                     

#  else 
#  if(method == "proportions" & !is.na(treat) & !is.na(control) & !is.na(n) & !is.na(power)){
#   # Maximum detectable difference. From Fleiss (1981).
#   z.beta <- qnorm(power, mean = 0, sd = 1) 
#   C <- (z.alpha + z.beta)^2
#   p <- ((treat * (1 - treat)) + (control * (1 - control)))
#   delta <- 2 / ((n - 2) - (C * p))
#   rval <- list(delta = delta)
#      }

#   else 
#   if(method == "proportions" & is.na(power)){
#   # Study power.  From Fleiss (1981).
#   delta <- abs(treat - control)
#   s1 <- delta^2 * (n - 2 - (2 / delta))
#   s2 <- ((treat * (1 - treat)) + (control * (1 - control))) 
#   z.beta <- sqrt(s1/s2) - z.alpha
#   power <- pnorm(z.beta, mean = 0, sd = 1)
#   rval <- list(power = power)
#      }

 else 
 if(method == "survival" & !is.na(treat) & !is.na(control) & is.na(n) & !is.na(power)){
  # Sample size.
  # From: Therneau TM and Grambsch PM 2000. Modelling Survival Data - Extending the Cox Model. Springer, London, p 61 - 65.
  z.beta <- qnorm(power, mean = 0, sd = 1)
  p <- r / (r + 1); q <- 1 - p
  # p <- 0.5; q <- 1 - p
  exp.beta <- log(treat) / log(control)
  n <- ((z.alpha + z.beta)^2) / (p * q * log(exp.beta)^2)
  
  # Account for the design effect:
  n <- n * design
  
  n.crude <- ceiling(n)
  n.treat <- ceiling(n / (r + 1)) * r
  n.control <- ceiling(n / (r + 1)) * 1
  n.total <- n.treat + n.control
  rval <- list(n.crude = n.crude, n.total = n.total, n.treat = n.treat, n.control = n.control)
  }

  else 
  if(method == "survival" & !is.na(treat) & !is.na(control) & !is.na(n) & is.na(power)){
  # Power.
  # From: Therneau TM and Grambsch PM 2000. Modelling Survival Data - Extending the Cox Model. Springer, London, p 61 - 65. 
  beta <- log(treat / control)
  p <- r / (r + 1); q <- 1 - p
  
  # Account for the design effect:
  n <- n / design
  
  z.beta <- sqrt(n * p * q * beta^2) - z.alpha
  power <- pnorm(z.beta, mean = 0, sd = 1)
  rval <- list(power = power)
     }
  
 else 
 if(method == "survival" & is.na(treat) & is.na(control) & !is.na(n) & !is.na(power)){
  # Maximum detectable difference.
  # From: Therneau TM and Grambsch PM 2000. Modelling Survival Data - Extending the Cox Model. Springer, London, p 61 - 65. 
  p <- r / (r + 1); q <- 1 - p
  z.beta <- qnorm(power, mean = 0, sd = 1) 
  
  # Account for the design effect:
  n <- n / design
  
  beta <- sqrt(((z.alpha + z.beta)^2) / (n * p * q))
  delta <- exp(beta)
  rval <- list(hazard = sort(c(delta, 1/delta)))
     }

 else 
 if(method == "cohort.count" & !is.na(treat) & !is.na(control) & is.na(n) & !is.na(power)){
  # Sample size estimate. From Woodward p 405:
  z.beta <- qnorm(power, mean = 0, sd = 1) 
  lambda <- treat / control
  pi <- control
  pc <- (pi * ((r * lambda) + 1)) / (r + 1)
  p1 <- (r + 1) / (r * (lambda - 1)^2 * pi^2)
  p2 <- z.alpha * sqrt((r + 1) * pc * (1 - pc))
  p3 <- z.beta * sqrt((lambda * pi * (1 - (lambda * pi))) + (r * pi * (1 - pi)))
  n <- p1 * (p2 + p3)^2
  
  # Account for the design effect:
  n <- n * design
  
  n.crude <- ceiling(n)
  n.treat <- ceiling(n / (r + 1)) * r
  n.control <- ceiling(n / (r + 1)) * 1
  n.total <- n.treat + n.control
  rval <- list(n.crude = n.crude, n.total = n.total, n.treat = n.treat, n.control = n.control)
  }

  else 
  if(method == "cohort.count" & !is.na(treat) & !is.na(control) & !is.na(n) & is.na(power)){
  # Study power. From Woodward p 409:
  lambda <- treat / control
  pi <- control
  pc <- (pi * ((r * lambda) + 1)) / (r + 1)

  # Account for the design effect:
  n <- n / design
  
  t1 <- ifelse(lambda >= 1, 
     (pi * (lambda - 1) * sqrt(n * r)),
     (pi * (1 - lambda) * sqrt(n * r)))
     
  t2 <- z.alpha * (r + 1) * sqrt(pc * (1 - pc))
  t3 <- (r + 1) * (lambda * pi * (1 - lambda * pi) + r * pi * (1 - pi))
  z.beta <- (t1 - t2) / sqrt(t3)
  power <- pnorm(z.beta, mean = 0, sd = 1)
  rval <- list(power = power)
     }

 else 
 if(method == "cohort.count" & is.na(treat) & !is.na(control) & !is.na(n) & !is.na(power)){
  # Risk ratio to be detected - requires a value for control. From Woodward p 409:
  z.beta <- qnorm(power, mean = 0, sd = 1) 
  pi <- control
  
  # Account for the design effect:
  n <- n / design
  
  Y <- r * n * pi^2
  Z <- (r + 1) * pi * (z.alpha + z.beta)^2
  a <- Y + (pi * Z)
  b <- (2 * Y) + Z
  c <- Y - (r * (1 - pi) * Z)
  lambda.pos <- (1 / (2 * a)) * (b + sqrt(b^2 - 4 * a * c))
  lambda.neg <- (1 / (2 * a)) * (b - sqrt(b^2 - 4 * a * c))
  rval <- list(lambda = sort(c(lambda.neg, lambda.pos)))
     }

 else 
 if(method == "case.control" & !is.na(treat) & !is.na(control) & is.na(n) & !is.na(power)){
  # Sample size. From Woodward p 412:
  z.beta <- qnorm(power, mean = 0, sd = 1) 
  lambda <- treat / control
  
  # For this function, 'sigma' is the proportion of study subjects exposed:
  Pc <- (sigma / (r + 1)) * ((r * lambda / (1 + (lambda - 1) * sigma)) + 1)
  T1 <- (r + 1) * (1 + (lambda - 1) * sigma)^2
  T2 <- r * sigma^2 * (sigma - 1)^2 * (lambda - 1)^2
  T3 <- z.alpha * sqrt((r + 1) * Pc * (1 - Pc))
  T4 <- lambda * sigma * (1 - sigma)
  T5 <- 1 + (lambda - 1) * sigma
  T6 <- T4 / (T5^2)
  T7 <- r * sigma * (1 - sigma)
  T8 <- z.beta * sqrt(T6 + T7)
  n <- (T1 / T2) * (T3 + T8)^2

  # P <- sigma
  # pc. <- (P / (r + 1)) * ((r * lambda) / (1 + ((lambda - 1) * P)) + 1)
  # p1 <- (r + 1) * (1 + (lambda - 1) * P)^2 / (r * P^2 * (P - 1)^2 * (lambda - 1)^2)
  # p2 <- z.alpha * sqrt((r + 1) * pc. * (1 - pc.))
  # p3 <- z.beta * sqrt(((lambda * P * (1 - P)) / ((1 + (lambda - 1) * P)^2)) + (r * P * (1 - P)))
  # n <- p1 * (p2 + p3)^2
  
  # Account for the design effect:
  n <- n * design
  
  n.crude <- ceiling(n)
  n.treat <- ceiling(n / (r + 1)) * r
  n.control <- ceiling(n / (r + 1)) * 1
  n.total <- n.treat + n.control
  rval <- list(n.crude = n.crude, n.total = n.total, n.treat = n.treat, n.control = n.control)
  }

 else 
 if(method == "case.control" & !is.na(treat) & !is.na(control) & !is.na(n) & is.na(power)){
  # Study power. From Woodward p 413:
  lambda <- treat / control
  # For this function, 'sd' is the proportion of study subjects exposed:
  P <- sigma
  # In this function "r" is input as the ratio of cases to controls. The formulae in Woodward assumes "r" is the ratio of controls to cases.
  r <- 1 /r
  
  # Account for the design effect:
  n <- n / design
  
  pc. <- (P / (r + 1)) * ((r * lambda) / (1 + ((lambda - 1) * P)) + 1)
  M <- abs(((lambda - 1) * (P - 1)) / (1 + (lambda - 1) * P))
  term.n1 <- (M * P * sqrt(n * r)) / sqrt(r + 1)
  term.n2 <- z.alpha * sqrt((r + 1) * pc. * (1 - pc.))
  term.d1 <- lambda * P * (1 - P) / (1 + (lambda - 1) * P)^2
  term.d2 <- r * P * (1 - P)
  z.beta <- (term.n1 - term.n2) / sqrt(term.d1 + term.d2)  
  power <- pnorm(z.beta, mean = 0, sd = 1)  
  rval <- list(power = power)
     }

  else 
  if(method == "case.control" & is.na(treat) & is.na(control) & !is.na(n) & !is.na(power)){
  # Risk ratio to be detected. From Woodward p 409:
  z.beta <- qnorm(power, mean = 0, sd = 1) 
  P <- sigma
  
  # Account for the design effect:
  n <- n / design
  
  a <- (r * P^2) - (n * r * P * (1 - P)) / ((z.alpha + z.beta)^2 * (r + 1))
  b <- 1 + (2 * r * P)
  lambda.pos <- 1 + ((-b + sqrt(b^2 - (4 * a * (r + 1)))) / (2 * a))
  lambda.neg <- 1 + ((-b - sqrt(b^2 - (4 * a * (r + 1)))) / (2 * a))
  rval <- list(lambda = sort(c(lambda.neg, lambda.pos)))
  }
rval
}
