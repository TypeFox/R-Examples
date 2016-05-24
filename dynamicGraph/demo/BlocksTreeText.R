
# Test with block recursive model:

# demo("startup.0", package = "dynamicGraph", verbose = FALSE)

source(paste(system.file(package = "dynamicGraph"), 
             "demo.source/startup.0.R", sep = "/"))

Block.tree <- list(Vertices = c("country"),
                   X = list(Vertices = c("sex", "race"),
                            A = list(Vertices = c("hair"),
                                     horizontal = TRUE, closed = TRUE,
                                     a = list(Vertices = c("eye"),
                                              horizontal = FALSE)),
                            B = list(Vertices = c("education")),
                            C = list(Vertices = c("age"))))

From <- c("country", "country",  "sex", "sex",      "race", "race")
To   <- c(    "sex",    "race", "hair", "eye", "education",  "age")

v <- unlist(Block.tree)
V.Names <- v[grep("Vertices", names(v))]
rm(v)
names(V.Names) <- NULL

V.Names <- c("alfa", V.Names, "z")

V.Types <- rep("Discrete", length(V.Names))

From <- match(From, V.Names)
To   <- match(To, V.Names)

control <- dg.control(UserMenus = Menus, updateEdgeLabels = FALSE,
                      useNamesForLabels = FALSE,
                      margin = 400, width = 600, height = 600, 
                      drawblocks = TRUE, overlaying = FALSE,
                      edgeColor = "green", vertexColor = "blue")

simpleGraph.Z.BTT <- new("dg.simple.graph", vertex.names = V.Names, 
                         types = V.Types, block.tree = Block.tree,
                         from = From, to = To,
                         # title = "<< BlocksTreeText.R >>",
                         texts = c("One\nlittle\ntext", 
                                   "An other\nlittle text"))

graph.Z.BTT <- simpleGraphToGraph(simpleGraph.Z.BTT, control = control)

Z.BTT <- dg(graph.Z.BTT, modelObject = Object, control = control,
            title = "BlocksTreeText.R")
