"quadline" <-
function (lm.obj, ...) 
{
x.range <- range(model.matrix(lm.obj)[,2])
x.df <- data.frame(seq(x.range[1],x.range[2],length=101))
names(x.df)<-names(coefficients(lm.obj))[2]
y.hat <- predict(lm.obj, newdata = x.df)
z <- spline(x.df[,1],y.hat)
lines(z$x, z$y, ...)
invisible()
}
