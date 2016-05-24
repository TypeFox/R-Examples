
# require(tcltk)

# demo("source")

# demo("usermenus", package = "dynamicGraph", verbose = FALSE)

# demo("defaultObjects", package = "dynamicGraph", verbose = FALSE)

source(paste(system.file(package = "dynamicGraph"), 
             "demo.source/defaultObjects.R", sep = "/"))

source(paste(system.file(package = "dynamicGraph"), 
             "demo.source/usermenus.R", sep = "/"))

# demo("colorlists", package = "dynamicGraph", verbose = FALSE)

# Object <- newYourModelObject("AnModelObject")

Object <- new("your.Model", name = "AnModelObject")

debug.strata   <- FALSE
debug.edges    <- FALSE
debug.position <- FALSE
debug.update   <- FALSE

V.Types <- c("Discrete", "Discrete", "Continuous", "Continuous")

V.Names <- c("Sex", "Age", "Eye", "Hair")
V.Names <- paste(V.Names, 1:4, sep ="/")

From <- c(1, 2, 3, 4)
To   <- c(2, 3, 4, 1)
