plot.mcsimex <-
function(
x, # Object of class mcsimex
xlab = expression((1 + lambda)), # name for the xlab
ylab = colnames(b[,-1]), # vector of labels for the yaxis
ask  = FALSE, # logical stopping after each plot
show = rep(TRUE, NCOL(b) - 1), # vector of logicals indicating what to plot
...) # optional parameters passed to par()
{
old.par <- par(no.readonly = TRUE)
on.exit(par(old.par))
par(...)
if (ask) par(ask = TRUE)
p.names <- names(coef(x))
b <- x$SIMEX.estimates
a <- seq(-1, max(b[,1]), by = 0.01)
d <- matrix(data=NA,nrow=length(a),ncol=NCOL(b)-1)
switch(x$fitting.method,
"quad" = d <- predict(x$extrapolation, newdata = data.frame(lambda = a)),
"line" = d <- predict(x$extrapolation, newdata = data.frame(lambda = a)),
"nonl" = for (i in 1:length(p.names))
d[, i] <- predict(x$extrapolation[[p.names[i]]],
newdata = data.frame(lambda = a)),
"log2" = for (i in 1:length(p.names))
d[, i] <- predict(x$extrapolation[[p.names[i]]],
newdata = data.frame(lambda = a)) -
((abs(apply(x$SIMEX.estimates[-1, -1], 2, min)) + 1) *
(apply(x$SIMEX.estimates[-1, -1], 2, min) <= 0))[i],
"logl" = d <- t(t(exp(predict(x$extrapolation,
newdata = data.frame(lambda = a)))) -
((abs(apply(x$SIMEX.estimates[-1, -1], 2, min)) + 1) *
(apply(x$SIMEX.estimates[-1, -1], 2, min) <= 0)))
)
for (i in 2:NCOL(b)) {
if (show[i - 1]) {
plot(b[, 1] + 1 , b[, i], xlab = xlab, ylab = ylab[i - 1], type = "n")
points(b[-1, 1] + 1, b[-1, i], pch = 19)
points(b[1, 1] + 1, b[1, i])
lines(a[a > 0] + 1, d[a > 0, (i - 1)])
lines(a[a < 0] + 1, d[a < 0, (i - 1)], lty = 2)
}
}
}

