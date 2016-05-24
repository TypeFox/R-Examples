
# Test, oriented (cyclic) edges, without causal structure:

# demo("startup", package = "dynamicGraph", verbose = FALSE)

source(paste(system.file(package = "dynamicGraph"), 
             "demo.source/startup.R", sep = "/"))

control <- dg.control(UserMenus = Menus, diagonal = TRUE, w = 4,
                      title = "<<Circle with extra vertices>>")

simpleGraph.Z.e <- new("dg.simple.graph", 
                       vertex.names = V.Names, types = V.Types, 
                       from = From, to = To, labels = V.Labels,
                       texts = c("-50-  0", "  0- 50", " 50-100", "100-150",
                                 "150-200", "200-250", "250-300", "300-350",
                                 "350-400", "400-450", "450-500"), 
                       extra.from = c(1:6, -(1:10)), 
                       extra.to = c(rep(-6, 6), -(2:11)))

# graph.Z.e <- as(simpleGraph.Z.e, "dg.graph")

graph.Z.e <- simpleGraphToGraph(simpleGraph.Z.e, diagonal = TRUE)

Z.e <- dg(graph.Z.e, modelObject = Object, control = control, 
          title = "Circle with extra vertices")

Z <- Z.e
