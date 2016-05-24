library(grid)
library(gridSVG)

# There are many filter effect primitives and many of them are
# quite complex so this file could end up having plenty more tests!

pdf(file = NULL)

# First, lets draw some text that we're then going to filter
grid.text("hello, world!", gp = gpar(fontsize = 96),
          name = "backtext")
# Draw a copy over the top with white text that will be left alone
grid.text("hello, world!", gp = gpar(fontsize = 96, col = "white"),
          name = "foretext")

# We want to create a filter that takes the text thicker, and then blurs it
f <- filterEffect(list(feMorphology(operator = "dilate",
                                    radius = unit(1, "mm")),
                       feGaussianBlur(sd = 1)))
# Apply the filter
grid.filter("backtext", f)

# Now lets export this 
grid.export("filter-test.svg")
dev.off()

