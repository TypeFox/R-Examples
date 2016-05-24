"plotNormXm" <-
function (x, im)
{
for (i in 1:im){
plot ( density(x[,i], na.rm=TRUE), col=4, main = paste("Density of ", colnames(x)[i]),
sub = "With Corresponding Normal" )
lines ( density(qnorm (ppoints(length(x[!(is.na(x))])),
mean(x[,i],na.rm=TRUE), sd(x[,i],na.rm=TRUE))), lty=2, col=2)
}
}

