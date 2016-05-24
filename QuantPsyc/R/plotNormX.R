"plotNormX" <-
function (x)
{
plot (density(x, na.rm=TRUE), col=4, main = paste("Density of ", deparse(substitute(x))),
sub = "With Corresponding Normal" )
lines ( density(qnorm (ppoints(length (x[!(is.na(x))])),
mean(x, na.rm=TRUE), sd(x, na.rm=TRUE))), lty=2, col=2)
}

