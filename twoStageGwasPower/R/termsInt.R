termsInt <-
function(aa, c.joint, pi.samples=0.4) {
  # prob z.joint > c.joint  givnen z.1 = aa  under null hypothesis
  mu.joint <- aa*sqrt(pi.samples)
  var.joint <- 1 - pi.samples
  sd.joint <- sqrt(var.joint)
  term.1 <- pnorm(c.joint, mean=mu.joint, sd=sd.joint, lower.tail=F)

  # prob z.joint < -c.joint  givnen z.1 = aa
  term.2 <- pnorm(-c.joint, mean=mu.joint, sd=sd.joint, lower.tail=T)
  
  result <- term.1 + term.2
  result2 <- list(term.1=term.1, term.2=term.2, result=result)
  result2
  }
