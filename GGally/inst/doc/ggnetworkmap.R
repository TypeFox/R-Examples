## ----echo=FALSE, message=FALSE-------------------------------------------
ignore <- suppressMessages(library(ggplot2))
ignore <- suppressMessages(library(maps))
ignore <- suppressMessages(library(network))
ignore <- suppressMessages(library(sna))
source("../R/utils.R")
source("../R/ggnetworkmap.R")
load("../data/twitter_spambots.rda")
knitr::opts_chunk$set(fig.width = 9, fig.height = 7, fig.retina = 1)

## ------------------------------------------------------------------------
airports <- read.csv("http://datasets.flowingdata.com/tuts/maparcs/airports.csv", header = TRUE)
rownames(airports) <- airports$iata

# select some random flights
set.seed(1234)
flights <- data.frame(
  origin = sample(airports[200:400, ]$iata, 200, replace = TRUE),
  destination = sample(airports[200:400, ]$iata, 200, replace = TRUE)
)

# convert to network
flights <- network(flights, directed = TRUE)

# add geographic coordinates
flights %v% "lat" <- airports[ network.vertex.names(flights), "lat" ]
flights %v% "lon" <- airports[ network.vertex.names(flights), "long" ]

# drop isolated airports
delete.vertices(flights, which(degree(flights) < 2))

# compute degree centrality
flights %v% "degree" <- degree(flights, gmode = "digraph")

# add random groups
flights %v% "mygroup" <- sample(letters[1:4], network.size(flights), replace = TRUE)

# create a map of the USA
usa <- ggplot(map_data("usa"), aes(x = long, y = lat)) +
  geom_polygon(aes(group = group), color = "grey65",
               fill = "#f9f9f9", size = 0.2)

delete.vertices(flights, which(flights %v% "lon" < min(usa$data$long)))
delete.vertices(flights, which(flights %v% "lon" > max(usa$data$long)))
delete.vertices(flights, which(flights %v% "lat" < min(usa$data$lat)))
delete.vertices(flights, which(flights %v% "lat" > max(usa$data$lat)))

# overlay network data to map
ggnetworkmap(usa, flights, size = 4, great.circles = TRUE,
             node.group = mygroup, segment.color = "steelblue",
             ring.group = degree, weight = degree)

## ---- eval=FALSE---------------------------------------------------------
#  data(twitter_spambots)

## ------------------------------------------------------------------------
# create a world map
world <- fortify(map("world", plot = FALSE, fill = TRUE))
world <- ggplot(world, aes(x = long, y = lat)) +
  geom_polygon(aes(group = group), color = "grey65",
               fill = "#f9f9f9", size = 0.2)

# view global structure
ggnetworkmap(world, twitter_spambots)

## ------------------------------------------------------------------------
ggnetworkmap(net = twitter_spambots, arrow.size = 0.5)

## ------------------------------------------------------------------------
# compute indegree and outdegree centrality
twitter_spambots %v% "indegree" <- degree(twitter_spambots, cmode = "indegree")
twitter_spambots %v% "outdegree" <- degree(twitter_spambots, cmode = "outdegree")

ggnetworkmap(net = twitter_spambots,
             arrow.size = 0.5,
             node.group = indegree,
             ring.group = outdegree, size = 4) +
  scale_fill_continuous("Indegree", high = "red", low = "yellow") +
  labs(color = "Outdegree")

## ------------------------------------------------------------------------
# show some vertex attributes associated with each account
ggnetworkmap(net = twitter_spambots,
             arrow.size = 0.5,
             node.group = followers,
             ring.group = friends,
             size = 4,
             weight = indegree,
             label.nodes = TRUE, vjust = -1.5) +
  scale_fill_continuous("Followers", high = "red", low = "yellow") +
  labs(color = "Friends") +
  scale_color_continuous(low = "lightgreen", high = "darkgreen")

