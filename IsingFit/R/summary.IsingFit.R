summary.IsingFit <-
function(object, ...)
{
  cat("\tNetwork Density:\t\t", round(mean(object$weiadj[upper.tri(object$weiadj)]!=0),2),"\n",
      "Gamma:\t\t\t",round(object$gamma,2),"\n",
      "Rule used:\t\t",ifelse(object$AND,"And-rule","Or-rule"),"\n",
      "Analysis took:\t\t",format(object$time,format="%s"),"\n"
      )
}
