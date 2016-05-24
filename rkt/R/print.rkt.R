print.rkt <-
function(x,...)
{
if (is.na(x$S))
{
cat("\nAnalysis not performed")
}
else
{
cat("\nStandard model")
cat("\nTau =",x$tau)
cat("\nScore = ",x$S)
cat("\nvar(Score) = ",x$varS)
cat("\n2-sided p-value = ",x$sl)
cat("\nTheil-Sen's (MK) or seasonal/regional Kendall (SKT/RKT) slope= ",x$B)
if (!is.na(x$varS.corrected))
{
cat("\n\nCorrection for inter-block covariance")
cat("\nvar(Score) = ",x$varS.corrected)
cat("\n2-sided p-value = ",x$sl.corrected)
}
if (!is.na(x$partial.varS))
{
cat("\n\nPartial model")
cat("\nPartial score = ",x$partial.S)
cat("\nvar(Partial score) = ",x$partial.varS)
cat("\n2-sided p-value = ",x$partial.sl)
if (!is.na(x$partial.varS.corrected))
{
cat("\n\nCorrection for inter-block covariance")
cat("\nvar(Score) = ",x$partial.varS.corrected)
cat("\n2-sided p-value = ",x$partial.sl.corrected)
}
}
}
cat("\n")
}
