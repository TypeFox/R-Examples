cec.plot.cost.function <- function(C, xlab="Iteration", ylab="Cost function", lwd=5, col="red", lwd.points=5, pch.points=19, col.points="black", mgp=c(1.5,0.5,0), ...)
{
    plot(x = 1:(length(C$cost) - 1), y = C$cost[2:(length(C$cost))], xlab=xlab, ylab=ylab, type="l", lwd=lwd, col=col, mgp=mgp, ...)
    points(x = 1:(length(C$cost) - 1), y = C$cost[2:(length(C$cost))], lwd=lwd.points, pch=pch.points, col=col.points)
    title("Cost function at each iteration")
}