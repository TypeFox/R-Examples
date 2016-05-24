`circ` <-
function()
{
# net(1)
C = RPMG::circle(1)
plot(C$x,C$y, type='n', axes=FALSE, asp=1, xlab="", ylab="")
 lines(c(C$x, C$x[1]),c(C$y,C$y[1]) , type='l')
}

