
# Test with block recursive model:

# demo("startup.0", package = "dynamicGraph", verbose = FALSE)

source(paste(system.file(package = "dynamicGraph"), 
             "demo.source/startup.0.R", sep = "/"))

Block.tree <- list(Vertices = c("country"), g = 30,
                   X = list(Vertices = c("sex", "race"), g = 30,
                            A = list(Vertices = c("hair", "eye"), g = 30,
                                     horizontal = FALSE),
                            B = list(Vertices = c("education"), g = 30),
                            C = list(Vertices = c("age"), g = 30)))

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
                      drawblocks = TRUE, overlaying = TRUE,
                      edgeColor = "green", vertexColor = "blue")

simpleGraph.Z.BTTG <- new("dg.simple.graph", vertex.names = V.Names, 
                          types = V.Types, block.tree = Block.tree,
                          from = From, to = To, 
                          # label = "<< BlocksTreeText.R >>",
                          texts = c("One\nlittle\ntext", 
                                    "An other\nlittle text"))

graph.Z.BTTG <- simpleGraphToGraph(simpleGraph.Z.BTTG, control = control)

Z.BTTG <- dg(graph.Z.BTTG, modelObject = Object, control = control, 
            title = "BlocksTreeTextG.R")
