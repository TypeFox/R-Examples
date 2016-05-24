twoStagePower <-
function(pi.samples, pi.markers, p0, p1, c1, c2, c.joint, c.singleStage, n.cases, n.controls) {
  # compute power given allele frequencies p0 and p1 in controls and cases
  mu1 <- (p1 - p0)/sqrt((p1*(1-p1)/n.cases + p0*(1-p0)/n.controls)/(2*pi.samples))
  mu.singleStage <- (p1 - p0)/sqrt((p1*(1-p1)/n.cases + p0*(1-p0)/n.controls)/2)
  sd1 <- 1
  p.stageOne <- pnorm(q=c1, mean=mu1, sd=sd1, lower.tail=F) +
        pnorm(q=-c1, mean=mu1, sd=sd1, lower.tail=T)
  mu2 <- (p1 - p0)/sqrt((p1*(1-p1)/n.cases + p0*(1-p0)/n.controls)/(2*(1 - pi.samples)))
  sd2 <- 1
  p.stageTwo <- pnorm(q=c2, mean=mu2, sd=sd2, lower.tail=F)*pnorm(q=c1, mean=mu1, sd=sd1, lower.tail=F)/p.stageOne +
        pnorm(q=-c2, mean=mu2, sd=sd2, lower.tail=T)*pnorm(q=-c1, mean=mu1, sd=sd1, lower.tail=T)/p.stageOne
  power.rep <- p.stageOne*p.stageTwo
  p.joint <- find.p.joint.alt(c.joint=c.joint, c1=c1, pi.samples=pi.samples, pi.markers=pi.markers, 
                           p0=p0, p1=p1, n.cases=n.cases, n.controls=n.controls) 
  power.joint <- p.stageOne*p.joint
  
  power.singleStage <- pnorm(q=c.singleStage, mean=mu.singleStage, sd=sd1, lower.tail=F) +
        pnorm(q=-c.singleStage, mean=mu.singleStage, sd=sd1, lower.tail=T)
  result <- list(p.stageOne=p.stageOne, 
                 power.rep=power.rep, power.joint=power.joint, power.singleStage=power.singleStage)
  result
  }
