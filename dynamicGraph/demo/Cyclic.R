
# Test, oriented (cyclic) edges, without causal structure:

# demo("startup", package = "dynamicGraph", verbose = FALSE)

source(paste(system.file(package = "dynamicGraph"), 
             "demo.source/startup.R", sep = "/"))

control <- dg.control(UserMenus = Menus, width = 200, height = 200, 
                      margin = 200, title = "<<Oriented edges>>")

simpleGraph.Z.o <- new("dg.simple.graph", 
                       vertex.names = V.Names, types = V.Types, 
                       from = From, to = To, oriented = TRUE)

graph.Z.o <- as(simpleGraph.Z.o, "dg.graph")

Z.o <- dg(graph.Z.o, modelObject = Object, control = control, 
          title = "Oriented edges")

