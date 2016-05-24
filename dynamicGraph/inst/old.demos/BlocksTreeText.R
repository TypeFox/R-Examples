
# Test with block recursive model:

# demo("startup.0", package = "dynamicGraph", verbose = FALSE)

source(paste(system.file(package = "dynamicGraph"), 
             "Demo/startup.0.R", sep = "/"))

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

V.Types <- rep("Discrete", length(V.Names))

From <- match(From, V.Names)
To   <- match(To, V.Names)

Z <- DynamicGraph(V.Names, V.Types, From, To, block.tree = Block.tree,
                  object = Object, margin = 400,
                  width = 600, height = 600, drawblocks = TRUE,
                  UserMenus = Menus, overlaying = FALSE,
                  debug.strata = debug.strata, debug.edges = debug.edges, 
                  debug.position = debug.position, debug.update = debug.update)
