## ---- eval=TRUE, warning=FALSE, message=FALSE----------------------------
library(alakazam)

# Load Change-O file
file <- system.file("extdata", "ExampleDb.gz", package="alakazam")
df <- readChangeoDb(file)

# Select clone
sub_df <- subset(df, CLONE == 164)

## ---- eval=TRUE----------------------------------------------------------
# This example data set does not have ragged ends
# Preprocess clone without ragged end masking (default)
clone <- makeChangeoClone(sub_df, text_fields=c("SAMPLE", "ISOTYPE"), 
                          num_fields="DUPCOUNT")

# Show combined annotations
clone@data[, c("SAMPLE", "ISOTYPE", "DUPCOUNT")]

## ---- eval=TRUE, warning=FALSE, message=FALSE----------------------------
library(igraph)
# Run PHYLIP and parse output
dnapars_exec <- "~/apps/phylip-3.69/dnapars"
graph <- buildPhylipLineage(clone, dnapars_exec, rm_temp=TRUE)

# Clone annotations
data.frame(CLONE=graph$clone,
           JUNCTION_LENGTH=graph$junc_len,
           V_GENE=graph$v_gene,
           J_GENE=graph$j_gene)

# Sequence annotations
data.frame(SEQUENCE_ID=V(graph)$name, 
           ISOTYPE=V(graph)$ISOTYPE,
           DUPCOUNT=V(graph)$DUPCOUNT)

## ---- eval=TRUE----------------------------------------------------------
# Plot graph with defaults
plot(graph)

## ---- eval=TRUE----------------------------------------------------------
# Modify graph and plot attributes
V(graph)$color <- "lightgrey"
V(graph)$color[V(graph)$name == "Germline"] <- "black"
V(graph)$color[grepl("Inferred", V(graph)$name)] <- "white"
V(graph)$label <- V(graph)$ISOTYPE
E(graph)$label <- ""

# Remove large default margins
par(mar=c(0, 0, 0, 0) + 0.1)
# Define a tree layout with the Germline at the top
ly <- layout_as_tree(graph, root="Germline", circular=F, flip.y=T)
# Plot graph
plot(graph, layout=ly, edge.arrow.mode=0, vertex.frame.color="black",
     vertex.label.color="black", vertex.size=50)
# Add legend
legend("topleft", c("Germline", "Inferred", "Sample"), 
       fill=c("black", "white", "grey80"), cex=0.75)

## ---- eval=TRUE, warning=FALSE, results="hide"---------------------------
library(dplyr)

# Preprocess clones
clones <- df %>%
    group_by(CLONE) %>%
    do(CHANGEO=makeChangeoClone(., text_fields=c("SAMPLE", "ISOTYPE"), 
                                num_fields="DUPCOUNT"))

# Build lineages
dnapars_exec <- "~/apps/phylip-3.69/dnapars"
graphs <- lapply(clones$CHANGEO, buildPhylipLineage, 
                 dnapars_exec=dnapars_exec, rm_temp=TRUE)

## ---- eval=TRUE----------------------------------------------------------
# Note, clones with only a single sequence will not be processed.
# A warning will be generated and NULL will be returned by buildPhylipLineage
# These entries may be removed for clarity
graphs[sapply(graphs, is.null)] <- NULL

# Leaving a subset of clones
nrow(clones)
length(graphs)

