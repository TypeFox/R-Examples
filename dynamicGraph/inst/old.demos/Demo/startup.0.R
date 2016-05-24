
# require(tcltk)

# demo("source")

# demo("usermenus", package = "dynamicGraph", verbose = FALSE)

# demo("defaultObjects", package = "dynamicGraph", verbose = FALSE)

source(paste(system.file(package = "dynamicGraph"), 
             "Demo/defaultObjects.R", sep = "/"))

source(paste(system.file(package = "dynamicGraph"), 
             "Demo/usermenus.R", sep = "/"))

# demo("colorlists", package = "dynamicGraph", verbose = FALSE)

Object <- newYourModelObject("AnModelObject")

debug.strata   <- TRUE
debug.edges    <- TRUE
debug.position <- TRUE
debug.update   <- TRUE

debug.strata   <- FALSE
debug.edges    <- FALSE
debug.position <- FALSE
debug.update   <- FALSE
