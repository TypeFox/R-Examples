library(grid)
library(gridSVG)

pdf(file = NULL)

# Create the definition of an opacity mask, in this case
# it will just be a circle with a gradient fill.
rg <- radialGradient(col = c("white", "black"))
m <- mask(gradientFillGrob(circleGrob(), gradient = rg))

# Register the opacity mask so that we can refer to it and apply it.
registerMask("circleMask", m)

# Creating a simple plot that will be masked to the circle
sp <- gTree(name = "simplePlot",
            children = gList(xaxisGrob(), yaxisGrob(),
                             rectGrob(gp = gpar(fill = "grey"))),
            vp = plotViewport())
grid.draw(sp)

# Now lets mask it
grid.mask("simplePlot", label = "circleMask")
# Alternatively we could also do this which avoids the need to call
# 'registerMask' explicitly
# grid.mask("simplePlot", m)

grid.export("mask-simpleplot.svg")
dev.off()


# Now lets recreate the previous example using 'pushMask'
pdf(file = NULL)
# Clear previous mask reference
gridSVG.newpage()


# Create the definition of an opacity mask, in this case
# it will just be a circle with a gradient fill.
rg <- radialGradient(col = c("white", "black"))
m <- mask(gradientFillGrob(circleGrob(), gradient = rg))

# Create a new masking context for this viewport
pushMask(m)

# Creating a simple plot that will be masked to our current masking context
sp <- gTree(name = "simplePlot",
            children = gList(xaxisGrob(), yaxisGrob(),
                             rectGrob(gp = gpar(fill = "grey"))),
            vp = plotViewport())
grid.draw(sp)

# End the masking context
popMask()

grid.export("pushmask-simpleplot.svg")
dev.off()


