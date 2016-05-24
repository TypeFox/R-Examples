anova.mean <- function(...)
  .Defunct("anovaMean", package="HH")

"anovaMean" <-
function(object, n, ybar, s, ..., ylabel="ylabel") {
  ## the object argument contains the levels of the factor
  df <- sum(n)-length(ybar)
  sigmahat <- (sum((n-1)*s^2) / df)^.5

  SS.group <- sum(n*(ybar - (ybar %*% n)/sum(n))^2)

  Df <- c(length(ybar)-1, sum(n)-length(ybar))
  Sum.of.Sq <- c(SS.group, df*sigmahat^2)
  Mean.Sq <- c(SS.group/(length(ybar)-1), sigmahat^2)
  F.Value <- c(Mean.Sq[1]/Mean.Sq[2], NA)
  Pr.F <- c(1-pf(F.Value[1], Df[1], Df[2]), NA)
  result <- list(Df=Df, "Sum of Sq"=Sum.of.Sq, "Mean Sq"=Mean.Sq,
                 "F value"=F.Value, "Pr(F)"=Pr.F)
  oldClass(result) <- c("anova","data.frame")
  attr(result,"row.names") <- format(c(ylabel,"Residuals"))
  attr(result, "heading") <-
    c("Analysis of Variance Table\n",
      paste("Response: ", ylabel, "\n", sep=""),
      "Terms added sequentially (first to last)")
  result
}
