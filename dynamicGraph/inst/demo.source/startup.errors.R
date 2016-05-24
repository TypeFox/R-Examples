
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

V.Types <- c("Discrete", "Ordinal", "Discrete",
             "Continuous", "Discrete", "Continuous")

V.Names <- c("Sex", "Age", "Eye", "FEV", "Hair", "Shosize")
V.Names <- paste(V.Names, 1:6, sep ="/")

From <- c(V.Names, "Gaf")
To   <- c(V.Names[-1], V.Names[1], "Gaf")
To   <- c(V.Names[-1], "Gryf",    "Gryf")
