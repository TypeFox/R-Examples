"eda.uni" <-
function(x,title="")
{
par (mfrow = c(2,2))
hist(x,main=title, xlab=deparse(substitute(x)))
plot(density(x, na.rm=TRUE),main="Smoothed Histogram")
qqnorm(x); qqline(x)
boxplot(x,horizontal=TRUE,main="BoxPlot")
}

