
# Test with 6 edges, no causal structure:

# demo("startup", package = "dynamicGraph", verbose = FALSE)

source(paste(system.file(package = "dynamicGraph"), 
             "Demo/startup.R", sep = "/"))

Z <- DynamicGraph(V.Names, V.Types, From, To, labels = V.Labels,
                  texts = c("-50-  0", "  0- 50", " 50-100", "100-150",
                            "150-200", "200-250", "250-300", "300-350",
                            "350-400", "400-450", "450-500"), diagonal = TRUE, 
                  extra.from = c(1:6, -(1:10)), 
                  extra.to = c(rep(-6, 6), -(2:11)),
                  object = Object, UserMenus = Menus, w = 4,
                  debug.strata = debug.strata, debug.edges = debug.edges, 
                  debug.position = debug.position, debug.update = debug.update)
