arcAngle <-
function(x, y, center)
{
begin <- atan2(y[1] - center[2], x[1] - center[1])
end <- atan2(y[2] - center[2], x[2] - center[1])
if (begin > end) begin <- begin - 2 * pi

c(begin=begin, end=end)
}
