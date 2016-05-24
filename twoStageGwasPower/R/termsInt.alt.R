termsInt.alt <-
function(aa, p0, p1, c.joint, pi.samples, n.cases, n.controls) {
  # prob z.joint > c.joint  givnen z.1 = aa  under alternative hypothesis
  mu1 <- (p1 - p0)/sqrt((p1*(1-p1) + p0*(1-p0))/(2*n.cases*pi.samples))
  #mu.jointA <- (p1 - p0)/sqrt((p1*(1-p1) + p0*(1-p0))/(2*n.cases)) + sqrt(pi.samples)*(aa - mu1)
  mu.joint <- sqrt(1-pi.samples)*(p1 - p0)/sqrt((p1*(1-p1)/n.cases + p0*(1-p0)/n.controls)/(2*(1 - pi.samples))) + sqrt(pi.samples)*aa
  var.joint <- 1 - pi.samples
  sd.joint <- sqrt(var.joint)
  
  term.1 <- pnorm(c.joint, mean=mu.joint, sd=sd.joint, lower.tail=F)

  # prob z.joint < -c.joint  givnen z.1 = aa
  term.2 <- pnorm(-c.joint, mean=mu.joint, sd=sd.joint, lower.tail=T)
  
  result <- term.1 + term.2
  result2 <- list(term.1=term.1, term.2=term.2, result=result)
  result2
  }
