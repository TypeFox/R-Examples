print.twoStageNull <-
function(x, ...) {
  cat("Threshold for one-stage design", x$c.singleStage,"\n")
  cat("Threshold for stage one", x$c1,"\n")
  cat("Threshold for replication design", x$c2,"\n")
  cat("Threshold for joint design", x$c.joint,"\n\n")
  }
