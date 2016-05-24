#' @demo Simple routing demo using the igraph-package.
#'
#' @details
#'   Simple routing using the igraph-package. In Munich, we want from
#'   the "Sendlinger Tor" to the north.
#'
#'   OSM planet dump is copied from
#'   http://www.dbs.ifi.lmu.de/cms/Project_PAROS.
#'
#'   Ensure that Osmosis is installed and in your PATH environment
#'   variable; see http://wiki.openstreetmap.org/wiki/Osmosis.

library("osmar")



### Download and extract data: #######################################

download.file("http://osmar.r-forge.r-project.org/muenchen.osm.gz",
              "muenchen.osm.gz")

system("gzip -d muenchen.osm.gz")



### Import subset based on bounding box: #############################

src <- osmsource_osmosis(file = "muenchen.osm",
                         osmosis = "osmosis")

muc_bbox <- center_bbox(11.575278, 48.137222, 3000, 3000)

muc <- get_osm(muc_bbox, src)
muc



### Reduce to highways: ##############################################

hways_muc <- subset(muc, way_ids = find(muc, way(tags(k == "highway"))))
hways <- find(hways_muc, way(tags(k == "name")))
hways <- find_down(muc, way(hways))
hways_muc <- subset(muc, ids = hways)
hways_muc


## Plot complete data and highways on top:
plot(muc)
plot_ways(hways_muc, col = "red", add = TRUE)


## Plot street map only:
plot_nodes(muc, pch = 19, cex = 0.1, col = "lightgray")
plot_ways(hways_muc, add = TRUE)
plot_nodes(hways_muc, add = TRUE, pch = 19, cex = 0.6)



### Define route start and end nodes: ################################

hway_start_node <- local({
  id <- find(muc, node(tags(v == "Sendlinger Tor")))[1]
  find_nearest_node(muc, id, way(tags(k == "highway")))
})
hway_start <- subset(muc, node(hway_start_node))

hway_end_node <- local({
  id <- find(muc, node(attrs(lon > 11.59 & lat > 48.150)))[1]
  find_nearest_node(muc, id, way(tags(k == "highway")))
})
hway_end <- subset(muc, node(hway_end_node))


## Add the route start and and nodes to the plot:
plot_nodes(hway_start, add = TRUE, col = "red", pch = 19, cex = 2)
plot_nodes(hway_end, add = TRUE, col = "blue", pch = 19, cex = 2)



### Create street graph and compute shortest route: ##################

gr_muc <- as_igraph(hways_muc)
summary(gr_muc)


### Compute shortest route:

route <- get.shortest.paths(gr_muc,
                            from = as.character(hway_start_node),
                            to = as.character(hway_end_node))[[1]]

route_nodes <- as.numeric(V(gr_muc)[route]$name)

route_ids <- find_up(hways_muc, node(route_nodes))
route_muc <- subset(hways_muc, ids = route_ids)
route_muc


## Add route to the plot:
plot_nodes(route_muc, add = TRUE, col = "green", pch = 19)
plot_ways(route_muc, add = TRUE, col = "green", lwd = 2)



### Compute route details: ###########################################

## Node and way IDs:
node_ids <- route_muc$nodes$attrs$id

way_ids <- local({
  w <- match(node_ids, route_muc$ways$refs$ref)
  route_muc$ways$refs$id[w]
})


## Way names:
way_names <- local({
  n <- subset(route_muc$ways$tags, k == "name")
  n[match(way_ids, n$id), "v"]
})


## Node coordinates, distances and bearings:
node_coords <- route_muc$nodes$attrs[, c("lon", "lat")]

node_dirs <- local({
  n <- nrow(node_coords)
  from <- 1:(n-1)
  to <- 2:n

  cbind(dist = c(0,
        distHaversine(node_coords[from, ],
                      node_coords[to, ])),
        bear = c(0,
        bearing(node_coords[from, ],
                node_coords[to, ])))
})


## Route details:
compass <- function(bearing) {
  dir <- function(x) {
    switch(as.character(x),
           "0" = "N",
           "1" = "NNE",
           "2" = "NE",
           "3" = "ENE",
           "4" = "E",
           "5" = "ESE",
           "6" = "SE",
           "7" = "SSE",
           "8" = "S",
           "9" = "SSW",
           "10" = "SW",
           "11" = "WSW",
           "12" = "W",
           "13" = "WNW",
           "14" = "NW",
           "15" = "NNW",
           "16" = "N")
  }
  sapply(round(bearing / 22.5), dir)
}

route_details <- data.frame(way_names, node_dirs)
route_details$cdist <- cumsum(route_details$dist)
route_details$dir <- compass(route_details$bear)

route_details
