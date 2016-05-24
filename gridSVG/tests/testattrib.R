

library(grid)
library(gridSVG)

# A very simple test
dev.new(width=6, height=6)
# Test script chunk
grid.script(file="test.script")
# Some default settings
pushViewport(viewport(gp=gpar(col="black", fill=NA)))
grid.circle(r=0.1, gp=gpar(fill="red"), name="circgrob")
# Test setting SVG attribute
grid.garnish("circgrob", onclick="circle_click(evt)")
popViewport()

grid.export()
dev.off()


# Single attribute value on single grob
grid.newpage()
grid.circle(r=.1, gp=gpar(fill="black"), name="c")
grid.garnish("c", onmousedown="alert('ouch')")
grid.export("testattrcircle.svg")

# Multiple attribute values on single grob
grid.newpage()
pushViewport(viewport())
grid.points(1:3/4, 1:3/4, pch=c(1, 10, 16), name="p")
grid.garnish("p", 
             onmousedown=c("alert('pch=1')",
               "alert('pch=10')",
               "alert('pch=16')"),
             group=FALSE)
grid.export("testattrpoints.svg")

# Multiple garnishes (one with single value, one with multiple values)
grid.newpage()
grid.circle(x=1:3/4, r=.1, gp=gpar(fill="black"), name="c")
grid.garnish("c",
             onmouseover=c("alert('c1')", "alert('c2')", "alert('c3')"),
             group=FALSE)
grid.garnish("c",
             onmousedown="alert('click me!')")
grid.export("testmultattr.svg")
             
# Sneak an SVG attribute through via gpar()
grid.newpage()
grid.text("test",
          gp=gpar("text-decoration"="line-through"),
          name="tt")
grid.export("testsvgattr.svg")
