Vectors3D <- function (vectors) 
{
    hW = 0.6
    hL = 0.1
    value = 5
    hWstart = hW
    hLstart = hL
    valAux = value
    while (value != 0) {
        open3d(windowRect = c(100, 100, 800, 800))
        bg3d("white")
        CMx = vectors[, 1]
        CMy = vectors[, 2]
        CMz = vectors[, 3]
        CRx = vectors[, 4]
        CRy = vectors[, 5]
        CRz = vectors[, 6]
        module = sqrt((CRx - CMx) * (CRx - CMx) + (CRy - CMy) * 
            (CRy - CMy) + (CRz - CMz) * (CRz - CMz))
        Cx <- c(CMx, CRx)
        Cy <- c(CMy, CRy)
        Cz <- c(CMz, CRz)
        Ctx <- aperm(array(Cx, dim = c(length(Cx)/2, 2)))
        Cty <- aperm(array(Cy, dim = c(length(Cy)/2, 2)))
        Ctz <- aperm(array(Cz, dim = c(length(Cz)/2, 2)))
        plot3d(Cx, Cy, Cz, type = "p", box = FALSE, col = "darkblue", 
            size = 2.5, aspect = c(1, 1, 1), axes = FALSE, xlab = "", 
            ylab = "", zlab = "")
        for (i in 1:length(CRx)) {
            Arrows3D(c(CMx[i], CMy[i], CMz[i]), c(CRx[i], CRy[i], 
                CRz[i]), colo = "red", hL = hL, hW = hW, 
                transparent = FALSE, plane = "YZ")
        }
        rgl.material(color = "black")
        decorate3d()
        grid3d(c("x", "y+", "z"))
        while (valAux == value) {
            valAux <- Slider(value)
        }
        rgl.close()
        val <- 5 - valAux
        value <- valAux
        if (val != 0) {
            hW <- hWstart
            hL <- hLstart
            if (val < 0) {
                hW = hW * (Mod(val) * 2)
                hL = hL * (Mod(val) * 2)
            }
            else {
                hW = hW/(Mod(val) * 2)
                hL = hL/(Mod(val) * 2)
            }
        }
        else {
            hW <- hWstart
            hL <- hLstart
        }
    }
}
