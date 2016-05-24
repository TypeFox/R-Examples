
# Test, oriented (cyclic) edges, without causal structure:

# demo("startup", package = "dynamicGraph", verbose = FALSE)

source(paste(system.file(package = "dynamicGraph"), 
             "Demo/startup.R", sep = "/"))

Z <- DynamicGraph(V.Names, V.Types, From, To, oriented = TRUE, 
                  object = Object, UserMenus = Menus, variableFrame = FALSE, 
                  debug.strata = debug.strata, debug.edges = debug.edges, 
                  debug.position = debug.position, debug.update = debug.update)
