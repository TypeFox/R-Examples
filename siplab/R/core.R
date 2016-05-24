core <-
function(plants, distance) {
# TRUE for plants at more than the given distance from the edge
    inside.owin(x = plants, w = erosion(as.owin(plants), distance))
}
