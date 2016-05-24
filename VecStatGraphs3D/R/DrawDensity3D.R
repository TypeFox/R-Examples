DrawDensity3D <- function (vectors, Div = 40, Layers = 3, DrawAxes = FALSE) 
{	
    open3d(windowRect = c(100, 100, 800, 800))
    bg3d("white")
    Cx = vectors[, 1]
    Cy = vectors[, 2]
    Cz = vectors[, 3]
    Cr <- kde3d(x = Cx, y = Cy, z = Cz, n = Div)
    th <- seq(min(Cr$d), max(Cr$d), len = Layers + 2)
    ramp <- colorRamp(c("white", "yellow", "red"))
    colo <- rgb(ramp(seq(0, 1, length = Layers)), maxColorValue = 255)
    al <- seq(0.1, 0.6, len = Layers)
    module = sqrt(Cx * Cx + Cy * Cy + Cz * Cz)
    spheres3d(0, 0, 0, radius = max(module), color = "black", 
        front = "line", back = "line", lwd = 1, smooth = TRUE, 
        lit = TRUE, line_antialias = FALSE, alpha = 0.2)
    x <- c(0, max(module), 0, 0)
    y <- c(0, 0, max(module), 0)
    z <- c(0, 0, 0, max(module))
    labels <- c("", "X", "Y", "Z")
    i <- c(1, 2, 1, 3, 1, 4)
    text3d(x, y, z, labels, adj = 0.8, cex = 1.5, font = 2, color = "black")
    segments3d(x[i], y[i], z[i], lwd = 3)
    rgl.points(x = Cx, y = Cy, z = Cz, size = 3, color = "black")
    contour3d(Cr$d, level = th[c(-1, -(Layers + 2))], x = Cr$x, 
        y = Cr$y, z = Cr$z, alpha = al, color = colo, add = TRUE, 
        engine = "rgl", fill = TRUE, smooth = 2, material = "shiny")
    if (DrawAxes == TRUE) {
        axes3d()
    }
}
