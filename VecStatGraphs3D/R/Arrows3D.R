Arrows3D <- function (a, b, colo="SkyBlue3", hL=0.3, hW=0.5, 
		width=1, transparent=FALSE, plane="XY")
{
    if (plane == "XY") {
        hN <- c(0, 0, 1)
    }
    if (plane == "YZ") {
        hN <- c(1, 0, 0)
    }
    if (plane == "XZ") {
        hN <- c(0, 1, 0)
    }
    aL = sqrt(sum((b - a) * (b - a)))
    aU = (b - a)/aL
    hP = c(aU[2] * hN[3] - aU[3] * hN[2], -1 * (aU[1] * hN[3] - 
        aU[3] * hN[1]), aU[1] * hN[2] - aU[2] * hN[1])
    hP = hP/sqrt(sum(hP * hP))
    quads3d(rbind(b, b - aU * hL + hP * hL/2 * 
        hW, b - aU * hL, b - aU * hL - 
        hP * hL/2 * hW), color = colo, lit = transparent)
    lines3d(rbind(a, b - aU * hL), color = colo, lwd = width)
}
