rgl.sphcirc=function (CrossEq = 0, PeakDec = 0, radius = 1, deg=TRUE, col = "black", ...) {
if(deg){CrossEq=CrossEq*pi/180;PeakDec=PeakDec*pi/180}
    rotdata = rotate3d(cbind(cos(seq(-pi, pi, len = 100)), 
        0, sin(seq(-pi, pi, len = 100))), x = 1, y = 0, 
        z = 0, angle = (pi/2 - PeakDec))*radius
    rotdata = rotate3d(rotdata, x = 0, y = 0, z = 1, angle = -CrossEq)
    lines3d(rotdata, aspect = TRUE, col = col, ...)
}

