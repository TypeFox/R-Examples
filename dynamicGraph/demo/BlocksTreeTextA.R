
# Test with block recursive model:

# demo("startup.0", package = "dynamicGraph", verbose = FALSE)

source(paste(system.file(package = "dynamicGraph"), 
             "demo.source/startup.0.R", sep = "/"))

Block.tree <- list(Vertices = c("country"),
                   X = list(Vertices = c("sex", "race"),
                            A = list(Vertices = c("hair", "eye"),
                                     horizontal = FALSE),
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
                      width = 600, height = 600, 
                      drawblocks = TRUE, overlaying = FALSE,
                      edgeColor = "green", vertexColor = "blue")

simpleGraph.Z.BTTA <- new("dg.simple.graph", vertex.names = V.Names, 
                          types = V.Types, block.tree = Block.tree,
                          from = From, to = To, 
                          # label = "<< BlocksTreeTextA.R >>",
                          texts = c("One\nlittle\ntext", 
                                    "An other\nlittle text"))

graph.Z.BTTA <- simpleGraphToGraph(simpleGraph.Z.BTTA, control = control)

Z.BTTA <- dg(graph.Z.BTTA, modelObject = Object, control = control,
            title = "BlocksTreeTextA.R")
