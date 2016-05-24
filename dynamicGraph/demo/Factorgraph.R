
# Test, oriented (cyclic) edges, without causal structure:

# demo("startup", package = "dynamicGraph", verbose = FALSE)

source(paste(system.file(package = "dynamicGraph"), 
             "demo.source/startup.R", sep = "/"))

Factors <- list(c(1, 2, 3, 4), c(3, 4, 5), c(4, 5, 6))

control <- dg.control(UserMenus = Menus, width = 200, height = 200, 
                      margin = 200, title = "<<Factorgraph>>")

simpleGraph.Z.f <- new("dg.simple.graph", 
                       vertex.names = V.Names, types = V.Types, 
                       labels = V.Labels, factors = Factors)

graph.Z.f <- as(simpleGraph.Z.f, "dg.graph")

Z.f <- dg(graph.Z.f, modelObject = Object, control = control, 
          title = "Factorgraph")
