library(misc3d)

local({
    data(teapot)

    haveRGL <- suppressWarnings(require(rgl,quietly=TRUE))

    ttri <- makeTriangles(teapot$vertices, teapot$edges,
                          color = "red", color2 = "green")
    edges <- teapot$edges

    ttriDull <- updateTriangles(ttri,material="dull")
    ttriShiny <- updateTriangles(ttri,material="shiny")
    ttriMetal <- updateTriangles(ttri,material="metal")

    ## teapots with varying materials
    drawScene(ttri,screen=list(y=-30,x=40), scale = FALSE)
    drawScene(ttriDull,screen=list(y=-30,x=40), scale = FALSE)
    drawScene(ttriShiny,screen=list(y=-30,x=40), scale = FALSE)
    drawScene(ttriMetal,screen=list(y=-30,x=40), scale = FALSE)
    drawScene(ttri,screen=list(y=-30,x=40),lighting=perspLighting,
              scale = FALSE)
    if (haveRGL) drawScene.rgl(ttri)

    ## teapots with varying colors
    drawScene(updateTriangles(ttriShiny,color2=grey.colors(ncol(edges))),
              screen=list(y=-30,x=40), scale = FALSE)
    drawScene(updateTriangles(ttriMetal, color2=heat.colors(ncol(edges))),
              screen=list(y=-30,x=40), scale = FALSE)
    drawScene(updateTriangles(ttriMetal,color2=heat.colors(ncol(edges))),
              screen=list(y=-30,x=40), scale = FALSE, engine="grid")

    ## two teapots side by side
    hc <- heat.colors(ncol(edges))
    drawScene(list(updateTriangles(ttri, color2 = hc),
                   translateTriangles(ttri,z=4)),
              screen=list(y=-30,x=40), scale = FALSE)
    drawScene(list(updateTriangles(ttri, color2 = hc),
                   translateTriangles(ttri,z=4)),
              screen=list(y=-30,x=40), scale = FALSE, engine="grid")
    drawScene(list(updateTriangles(ttriShiny,color2=hc),
                   translateTriangles(ttriMetal,z=4)),
              screen=list(y=-30,x=40), scale = FALSE, engine="grid")
    if (haveRGL)
        drawScene.rgl(list(updateTriangles(ttri, color=hc),
                           translateTriangles(ttri,z=4)))

    ## nested teapots
    drawScene(list(updateTriangles(ttri,color="blue", fill=FALSE,
                                   col.mesh="blue"),
                   scaleTriangles(updateTriangles(ttriMetal, color2="red"),
                                  0.6)),
              screen=list(y=-30,x=20,y=-140), scale = FALSE)
    if (haveRGL)
        drawScene.rgl(list(updateTriangles(ttri, alpha = 0.5, color="blue"),
                           scaleTriangles(ttriMetal, 0.6)))

    ## teapot with smoothing (Phong shading)
    drawScene(updateTriangles(ttriMetal, color2 = hc, smooth = 1),
              screen=list(y=-30,x=40), scale = FALSE)
    drawScene(updateTriangles(ttriMetal, color2 = hc, smooth = 2),
              screen=list(y=-30,x=40), scale = FALSE)
    drawScene(updateTriangles(ttriMetal, color2 = hc, smooth = 3),
              screen=list(y=-30,x=40), scale = FALSE)
})
