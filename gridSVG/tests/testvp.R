
# Test that 'vp' slots [viewports (vpTrees, vpLists, vpStacks) and vpPaths]
# are being recorded correctly

library(grid)
library(gridSVG)

grid.newpage()
vp <- viewport(width=.5, height=.5, gp=gpar(fill="grey"))
grid.rect(vp=vp)
grid.export("grob-viewport.svg")

grid.newpage()
vp <- vpStack(viewport(width=.5, height=.5, gp=gpar(fill="grey")),
              viewport(gp=gpar(col="red")))
grid.rect(vp=vp)
grid.export("grob-vpStack.svg")

grid.newpage()
vp <- vpList(viewport(width=.5, height=.5, gp=gpar(fill="grey")),
             viewport(gp=gpar(col="red")))
grid.rect(vp=vp)
grid.export("grob-vpList.svg")

grid.newpage()
vp <- vpTree(viewport(width=.5, height=.5),
             vpList(viewport(width=.5, height=.5, gp=gpar(fill="grey")),
                    viewport(gp=gpar(col="red"))))
grid.rect(vp=vp)
grid.export("grob-vpTree.svg")

grid.newpage()
vp <- vpTree(viewport(width=.5, height=.5, name="p"),
             vpList(viewport(width=.5, height=.5, gp=gpar(fill="grey"),
                             name="c1"),
                    viewport(gp=gpar(col="red"),
                             name="c2")))
pushViewport(vp)
upViewport(0)
grid.rect(vp="p::c1")
grid.export("grob-vpPath.svg")

grid.newpage()
vp <- viewport(width=.5, height=.5, gp=gpar(fill="grey"))
grid.draw(gTree(children=gList(rectGrob()), vp=vp))
grid.export("gTree-viewport.svg")

# Pathological
grid.newpage()
vp <- vpTree(viewport(width=.5, height=.5, name="p"),
             vpList(viewport(width=.5, height=.5, gp=gpar(fill="grey"),
                             name="c1"),
                    viewport(gp=gpar(col="red"),
                             name="c2")))
pushViewport(vp)
upViewport(0)
grid.draw(gTree(childrenvp=vp,
                children=gList(rectGrob(vp="p::c1")),
                vp="p::c1"))
grid.export("gTree-vpPath.svg")
