cor.comb <-
function(y1, y2, y3, y4, x1, x2, x3, x4, sims=1000, hist=FALSE, rnd=5, CI=TRUE, simple=FALSE, seed=2) {
  sim.dist <- rep(0, sims)
  x.set <- data.frame(x1, x2, x3, x4)
  y.set <- data.frame(y1, y2, y3, y4)
  n1 <- valid.pairs(y1, x1)$Valid
  n2 <- valid.pairs(y2, x2)$Valid
  n3 <- valid.pairs(y3, x3)$Valid
  n4 <- valid.pairs(y4, x4)$Valid
  Zr1 <- fisherz(cor(y.set$y1,x.set$x1, use="pair"))
  Zr2 <- fisherz(cor(y.set$y2,x.set$x2, use="pair"))
  Zr3 <- fisherz(cor(y.set$y3,x.set$x3, use="pair"))
  Zr4 <- fisherz(cor(y.set$y4,x.set$x4, use="pair"))
  Comb <- fisherz2r((Zr1 + Zr2 + Zr3 + Zr4) / 4)
  WgtES <- fisherz2r(sum((n1-3)*Zr1, (n2-3)*Zr2, (n3-3)*Zr3, (n4-3)*Zr4)/sum(n1-3,n2-3,n3-3,n4-3))

  if (simple==T) {return(round(Comb,rnd))}

  if(seed!=F) {set.seed(seed)}
  for (i in 1:sims) {
    order1 <- sample(nrow(x.set), nrow(x.set), replace=FALSE)   #Generate a sample of random orders
    psuedo.x <- x.set[order(order1),]
    sim1 <- fisherz(cor(y.set$y1,psuedo.x$x1, use="pair"))
    sim2 <- fisherz(cor(y.set$y2,psuedo.x$x2, use="pair"))
    sim3 <- fisherz(cor(y.set$y3,psuedo.x$x3, use="pair"))
    sim4 <- fisherz(cor(y.set$y4,psuedo.x$x4, use="pair"))
    SimAvg <- fisherz2r((sim1 + sim2 + sim3 + sim4) / 4)
    sim.dist[i] <- SimAvg
    }

    if (hist==TRUE) {
    hist(sim.dist)
    }

  if (Comb > median(sim.dist)) {randp <- sum(sim.dist > Comb)*2 / sims}
  if (Comb < median(sim.dist)) {randp <- sum(sim.dist < Comb)*2 / sims}
  if (Comb == median(sim.dist)) {randp <- 1.0}
  CI5 <- quantile(sim.dist, .025) + Comb
  CI95 <- quantile(sim.dist, .975) + Comb
  normp <- 2*pnorm(-abs((Comb - mean(sim.dist)) / (sd(sim.dist))))
    z1 <- qnorm(cor.test(y.set$y1,x.set$x1, use="pair")$p.value / 2)
    z2 <- qnorm(cor.test(y.set$y2,x.set$x2, use="pair")$p.value / 2)
    z3 <- qnorm(cor.test(y.set$y3,x.set$x3, use="pair")$p.value / 2)
    z4 <- qnorm(cor.test(y.set$y4,x.set$x4, use="pair")$p.value / 2)
    ifelse(Zr1 < 0, z1 <- -1*z1, ifelse(Zr1==0, z1 <- 0, z1 <- z1))
    ifelse(Zr2 < 0, z2 <- -1*z2, ifelse(Zr2==0, z2 <- 0, z2 <- z2))
    ifelse(Zr3 < 0, z3 <- -1*z3, ifelse(Zr3==0, z3 <- 0, z3 <- z3))
    ifelse(Zr4 < 0, z4 <- -1*z4, ifelse(Zr4==0, z4 <- 0, z4 <- z4))
  stoufp <- 2*pnorm(-abs((z1 + z2 + z3 + z4) / sqrt(4)))

  if (CI==F) {
  out <- round(cbind(Comb, WgtES, normp, randp, stoufp),rnd)
  colnames(out) <- c("UnWgt r", "Wgt r", "Normal p", "Rand p", "Stouffer p")
  rownames(out) <- c("Results")
  return(out)
  }

  else {
  out <- round(cbind(Comb, WgtES, normp, randp, stoufp, CI5, CI95),rnd)
  colnames(out) <- c("UnWgt r", "Wgt r", "Normal p", "Rand p", "Stouffer p", "LL", "UL")
  rownames(out) <- c("Results")
  return(out)
  }
}
