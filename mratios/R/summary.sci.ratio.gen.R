"summary.sci.ratio.gen" <-
function(object, digits=4,...)
{

cat("The following linear model has been fitted in lm: ","\n")

print(summary(object$fit), digits=digits, ...)

summary.sci.ratio(object, digits=digits, ...)


}

