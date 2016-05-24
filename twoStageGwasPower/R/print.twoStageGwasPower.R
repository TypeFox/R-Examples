print.twoStageGwasPower <-
function(x, ...) {
  cat("Power for one-stage design",x$power.singleStage,"\n")
  cat("Power for replication analysis", x$power.rep,"\n")
  cat("Power for joint analysis",x$power.joint,"\n\n")
  cat("Threshold for one-stage design", x$c.singleStage,"\n")
  cat("Threshold for stage one", x$c1,"\n")
  cat("Threshold for replication design", x$c2,"\n")
  cat("Threshold for joint design", x$c.joint,"\n\n")
  cat("Penetrance for GG", x$penetrance.GG,"\n")
  cat("Penetrance for Gg", x$penetrance.Gg,"\n")
  cat("Penetrance for gg", x$penetrance.gg,"\n\n")
  cat("Disease allele frequencies in cases", x$p1,"\n")
  cat("Disease allele frequencies in controls", x$p0,"\n")
  cat("Probability that associated markers will be followed up in stage two", x$p.stageOne, "\n\n") 
  cat("Reduction in genotyping from two-stage design", x$savings, "\n")
  }
