## ----------------------------------------------------------------------
##                   setSingleEnum
## ----------------------------------------------------------------------
ShapeEnum.gen <- setSingleEnum("Shape",
                                levels = c("circle", "line", "rectangle"))

obj <- new("ShapeSingleEnum", "circle")
obj
obj <- "triangle" # doesn't check, because it's not signal field.
obj # it's not SingleEnum object anymore, be careful.
class(obj) # just character

## only set it as properties, allow you to assign the value and
## validate it.
par.gen <- setRefClass("Graph",
                       properties(fields = list(shape = "ShapeSingleEnum"),
                                  prototype = list(shape = new("ShapeSingleEnum",
                                                     "circle"))))
pars <- par.gen$new()
pars$shape
pars$shape <- "line"
pars$shape
class(pars$shape)# still a SingleEnum
## ----------------------------------------------------------------------
##                   setMultipleEnum
## ----------------------------------------------------------------------
ShapeEnum.gen <- setMultipleEnum("Shape",
                                levels = c("circle", "line", "rectangle"))

par.gen <- setRefClass("Graph",
                       properties(list(shape = "ShapeMultipleEnum")))
## we can initialize in this way too
pars <- par.gen$new(shape = new("ShapeMultipleEnum", c("circle", "line")))
pars$shape
pars$shape <- c("line", "rectangle")
pars$shape
class(pars$shape)# still a MultipleEnum

## Color Single Enum
bgColorSingleEnum.gen <- setColorEnum("bgColor", levels = c("black", "white", "gray"))
obj <- new("bgColorSingleEnum", "white")
## Glyph Single Enum
PointSizeSingleEnum.gen <- setGlyphEnum("PointSize", levels = c("1", "2", "5", "10"), contains = "GlyphEnum")
obj <- new("PointSizeSingleEnum", "1")
obj
