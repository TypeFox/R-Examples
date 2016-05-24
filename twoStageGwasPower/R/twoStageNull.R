twoStageNull <-
function(pi.samples, pi.markers, alpha.marker) {
  # calculate thresholds
  # alpha.marker = alpha.genome/num.markers
  # typically, alpha.genome = 0.05
  c1 <- -qnorm(pi.markers/2)
  c2 <- -qnorm(alpha.marker/pi.markers)
  c.singleStage <- -qnorm(alpha.marker/2)
  result.opt <- optimize(f=find.p.joint.alpha, lower=2, upper=8, c1=c1, 
              pi.samples=pi.samples, pi.markers=pi.markers, alpha.marker=alpha.marker )
  c.joint <- result.opt$minimum
  result <- list( c1=c1, c2=c2, c.joint=c.joint, c.singleStage=c.singleStage)
  class(result) <- "twoStageNull"
  result
  }
