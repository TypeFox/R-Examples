library(grid)
library(gridSVG)

pdf(file = NULL)

# Create the definition of a clipping path, in this case
# it will just be a circle.
cp <- clipPath(circleGrob(r = 0.3))

# Register the clipping path so that we can refer to it and apply it.
registerClipPath("circleClip", cp)

# Creating a simple plot that will be clipped to the circle
sp <- gTree(name = "simplePlot",
            children = gList(xaxisGrob(), yaxisGrob(),
                             rectGrob(gp = gpar(fill = "grey"))),
            vp = plotViewport())
grid.draw(sp)

# Now lets clip it
grid.clipPath("simplePlot", label = "circleClip")
# Alternatively we could also do this which avoids the need to call
# 'registerClipPath' explicitly
# grid.clipPath("simplePlot", cp)

# All that remains is a grey circle
grid.export("clippath-simpleplot.svg")
dev.off()


# Now lets recreate the previous example using 'pushClipPath'
pdf(file = NULL)
# Clear previous clipping path reference
gridSVG.newpage()

pushViewport(plotViewport())

# Create the definition of a clipping path, in this case
# it will just be a circle.
cp <- clipPath(circleGrob(r = 0.3))

# Create a new clipping context for this viewport
pushClipPath(cp)

# Creating a simple plot that will be clipped to our current clipping context
sp <- gTree(name = "simplePlot",
            children = gList(xaxisGrob(), yaxisGrob(),
                             rectGrob(gp = gpar(fill = "grey"))))
grid.draw(sp)

# End the clipping context
popClipPath()
popViewport()

# Again we are just left with a grey circle
grid.export("pushclippath-simpleplot.svg")
dev.off()

