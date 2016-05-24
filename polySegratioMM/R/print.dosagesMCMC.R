`print.dosagesMCMC` <-
function(x, ..., index.sample=20)
{
  cat("Dosages for chain:", x$chain,"\n")
  cat("Thresholds set at:\n")
  print(x$thresholds)

  if (length(x$index.sample)==1) {
    index <- sort(sample(1:dim(x$p.dosage)[1],x$index.sample))
    cat("A random sample of posterior probabilities and classifications\n")
  } else {
    cat("Posterior probabilities and classifications for specified markers\n")
    index <- index.sample
  }      
  print(cbind(x$p.dosage[index,],x$dosage[index,],
              "maxPostP"=x$max.post.dosage[index]), na.print=".")
  ##,
  ##      quote=FALSE, na.print=".", zero.print = ".", justify="right")

  cat("\nMaximum posterior probabilities for",dim(x$p.dosage)[1],"markers\n")
  print(summary(x$max.post))
  cat("\nProportion of genes classified using maximum posterior",
      "probability\n")
  print(colMeans(x$p.dosage==x$max.post))
  cat("Total proportion of markers classified:",
      sum(colMeans(x$p.dosage==x$max.post)),"\n")

  cat("Call:\n")
  print(x$call)

}

