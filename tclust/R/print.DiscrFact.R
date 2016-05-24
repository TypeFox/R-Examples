print.DiscrFact<-
function (x, ...)
{
  cat ("Mean overall discriminant factor:", mean (x$assignfact), "\n")
#  cat ("\nFurther information on this DiscrFact - object goes here\n")

  cat ("Mean discriminant factor per cluster:\n")
  print (x$mean.DiscrFact)

  idx = x$assignfact > x$threshold

  if (!sum (idx))
  {
    cat ("No decision is considered as doubtful\n")
    return
  }

  cat (sum (idx), "decisions are considered as doubtful\n")
}
