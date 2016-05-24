## ----, setup, include=FALSE, echo=FALSE----------------------------------
library(intergraph)
library(knitr)
opts_knit$set( upload.fun = image_uri )
opts_chunk$set( highlight=TRUE )
options(markdown.HTML.options= unique(c(getOption("markdown.HTML.options"), "toc")))
set.seed(123)

## ----,packages-----------------------------------------------------------
library(intergraph)
library(network)
library(igraph)

## ----, summarize-igraph--------------------------------------------------
summary(exIgraph)
summary(exIgraph2)

## ----,summarize-network--------------------------------------------------
exNetwork
exNetwork2

## ----,network2igraph-----------------------------------------------------
# check class of 'exNetwork'
class(exNetwork)
# convert to 'igraph'
g <- asIgraph(exNetwork)
# check class of the result
class(g)

## ------------------------------------------------------------------------
el.g <- get.edgelist(g)
el.n <- as.matrix(exNetwork, "edgelist")
identical( as.numeric(el.g), as.numeric(el.n))

## ----,igraph2network-----------------------------------------------------
net <- asNetwork(exIgraph)

## ------------------------------------------------------------------------
el.g2 <- get.edgelist(exIgraph)
el.n2 <- as.matrix(net, "edgelist")
identical( as.numeric(el.g2), as.numeric(el.n2))

## ----attrmap-defaults----------------------------------------------------
attrmap()

## ----attrmap-example-rules-----------------------------------------------
new_rule <- data.frame(type="vertex", fromcls="network", fromattr="na",
                       tocls="igraph", toattr=NA,
                       stringsAsFactors=FALSE)
# combine with the default rules
rules <- rbind( attrmap(), new_rule )
rules

## ----attrmap-example-----------------------------------------------------
(ig1 <- asIgraph(exNetwork))
(ig2 <- asIgraph(exNetwork, amap=rules))

# check if "na" was dropped
"na" %in% igraph::list.vertex.attributes(ig1)
"na" %in% igraph::list.vertex.attributes(ig2)

## ----asDF----------------------------------------------------------------
l <- asDF(exIgraph)
str(l)

## ----show-edgedb---------------------------------------------------------
l$edges

## ----show-vertexdb-------------------------------------------------------
l$vertexes

## ----fromdf--------------------------------------------------------------
z <- asNetwork(l$edges, directed=TRUE, l$vertexes)
z

## ----showdata-code,eval=FALSE--------------------------------------------
#  layout(matrix(1:4, 2, 2, byrow=TRUE))
#  op <- par(mar=c(1,1,2,1))
#  # compute layout
#  coords <- layout.fruchterman.reingold(exIgraph)
#  plot(exIgraph, main="exIgraph", layout=coords)
#  plot(exIgraph2, main="exIgraph2", layout=coords)
#  plot(exNetwork, main="exNetwork", displaylabels=TRUE, coord=coords)
#  plot(exNetwork2, main="exNetwork2", displaylabels=TRUE, coord=coords)
#  par(op)

## ----showdata-pic, ref.label="showdata-code",echo=FALSE,fig.height=10,fig.width=10----
layout(matrix(1:4, 2, 2, byrow=TRUE))
op <- par(mar=c(1,1,2,1))
# compute layout
coords <- layout.fruchterman.reingold(exIgraph)
plot(exIgraph, main="exIgraph", layout=coords)
plot(exIgraph2, main="exIgraph2", layout=coords)
plot(exNetwork, main="exNetwork", displaylabels=TRUE, coord=coords)
plot(exNetwork2, main="exNetwork2", displaylabels=TRUE, coord=coords)
par(op)

## ----, session_info------------------------------------------------------
sessionInfo()

