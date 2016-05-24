createCircle <-
function(center, r, n=40, begin=0, end=2*pi)
{
theta  = seq(begin, end, (end - begin)/n)
x <- cos(theta) * r + center[1]
y <- sin(theta) * r + center[2]

cbind(x,y)
}
