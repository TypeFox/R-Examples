drawArc <-
function(x, y, center, r, ...) 
{
angles <- arcAngle(x, y, center)
createCircle(center, r, begin=angles[1], end=angles[2], ...)
}
