# icon = matrix(ncol=2, byrow=T, c(-4.002, -1.055,
# -5.514, -2.879,
# -10.733, -0.866,
# -10.414, -0.034,
# -11.23, -0.368,
# -11.593, 0.403,
# -11.6, -0.316,
# -12.126, 0.448,
# -12.801, 0.485,
# -16.93, -2.269,
# -18.501, -10.804,
# -27.307, -15.841,
# -27.025, -16.694,
# -31.398, -17.556,
# -36.054, -17.201,
# -36.647, -20.228,
# -39.627, -24.427,
# -40.798, -25.91,
# -43.037, -34.095,
# -49.664, -37.22,
# -55.416, -37.028,
# -60.787, -40.4,
# -66.432, -39.55,
# -67.739, -34.556,
# -64.846, -24.931,
# -66.141, -16.343,
# -62.978, -16.472,
# -59.509, -13.296,
# -61.866, -10.283,
# -64.234, -2.822,
# -65.498, 4.661,
# -64.607, 8.948,
# -62.429, 14.761,
# -63.997, 20.149,
# -63.762, 22.944,
# -65.335, 25.571,
# -57.167, 28.401,
# -52.631, 30.339,
# -49.399, 28.721,
# -46.449, 30.495,
# -41.423, 30.799,
# -34.196, 28.625,
# -29.8, 25.809,
# -26.546, 21.443,
# -25.486, 22.095,
# -24.797, 19.431,
# -20.998, 16.087,
# -14.315, 14.161,
# -6.45, 5.52))
# icon.x = icon[,1]
# icon.y = icon[,2]
# lines(icon.x/120 +1.7, icon.y/120 + 1.7)
#  grid.path - nothing shows up.

## currently, my php files must be in /Library/WebServer/Documents
##  ~Roger/Sites is not working.
## php is working: see /Library/WebServer/Documents/SVGconvert.php
## which conputed these numbers.
## Code came from  http://www.stoimen.com/blog/2011/02/11/from-svg-to-geo-coordinates-a-complete-guide/

drawIcon = function(x=0, y=0, icon, pch=1, closepath=T, color='black') {
  require(testthat)
  icon.x = icon[,1]
  icon.y = icon[,2]
  if(closepath) {
    icon.x = c(icon.x, icon.x[1]); icon.y = c(icon.y, icon.y[1])
  }
  icon.x = (icon.x - mean(icon.x)) * pch + x
  icon.y = (icon.y - mean(icon.y)) * pch + y
  lines(icon.x, icon.y, col=color)
}
