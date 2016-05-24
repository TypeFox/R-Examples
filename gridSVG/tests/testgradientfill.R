library(grid)
library(gridSVG)

pdf(file = NULL)

# Define a linear gradient
lg <- linearGradient(col = c("blue", "red"))
# Register it 
registerGradientFill("lingrad", lg)

# Do the same thing for radial gradient but also set the focus for the
# radial fill to be off-centre
rg <- radialGradient(fx = 0.3, fy = 0.3,
                     col = c("white", "black"), stops = c(0, 2))
registerGradientFill("radgrad", rg)

# Create rects that we are going to be gradient filling
grid.rect(x = 0.2, width = 0.2, height = 0.2, name = "linearfill")
grid.rect(x = 0.8, width = 0.2, height = 0.2, name = "radialfill")

# Now apply the gradients
grid.gradientFill("linearfill", label = "lingrad", alpha = 0.7)
grid.gradientFill("radialfill", label = "radgrad")

grid.export("gradient-test.svg")
dev.off()
