
# Test of "link" and "add":

# demo("startup", package = "dynamicGraph", verbose = FALSE)

source(paste(system.file(package = "dynamicGraph"), 
             "demo.source/startup.R", sep = "/"))

control <- dg.control(UserMenus = Menus, returnLink = TRUE, 
                      width = 150, height = 200, margin = 200)

simpleGraph.Z.1 <- new("dg.simple.graph", 
                       vertex.names = V.Names[1:3], types = V.Types[1:3], 
                       from = c(1, 2, 3), to = c(2, 3, 1))

graph.Z.1 <- as(simpleGraph.Z.1, "dg.graph")

Z <- dg(graph.Z.1, modelObject = Object, control = control, title = "Z")


simpleGraph.WV.1 <- new("dg.simple.graph", from = c(2, 3), to = c(3, 1))

graph.WV.1 <- simpleGraphToGraph(simpleGraph.WV.1,
                               vertexList = graph.Z.1@vertexList,
                               blockList = graph.Z.1@blockList)

W <- addModel(graph.WV.1, 
              frameModels = Z, control = control, title = "W")

V <- addView(graph.WV.1, 
             frameModels = Z, modelIndex = 1, control = control,
             viewType = "Factor", title = "W")


simpleGraph.U.1 <- new("dg.simple.graph", from = 1, to = 2)

graph.U.1 <- simpleGraphToGraph(simpleGraph.U.1,
                               Vertices = graph.Z.1@vertexList,
                               BlockList = graph.Z.1@blockList)

U <- replaceModel(graph.U.1, 
                  frameModels = Z, modelIndex = 1, graphIndex = 1, 
                  control = control, title = "U")
