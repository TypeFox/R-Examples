library(spacom)


## create dummy locations and deduce a weight matrix
points <- list()
points[['a']] <- c(0,4)
points[['b']] <- c(0,0)
points[['c']] <- c(7,4)
points[['d']] <- c(4,1)
nb.area <- length(points)
dist <- matrix(ncol=nb.area, nrow=nb.area)
distance <- function(x,y) {
    d <- x-y;
    return(sqrt(d[1]^2+d[2]^2))
}

for (i in 1:nb.area) {
    for (j in 1:nb.area) {
        dist[i, j] <- distance(points[[i]], points[[j]])
    }
}
geow <- WeightMatrix(dist, bandwidth = 4)

## create dummy precise data
area.name <- names(points)
colnames(geow) <- area.name
prec.1 <- runif(nb.area, 1, 2)
prec.1 <- data.frame(area.name, prec.1)
names(prec.1) <- c("area.name", "prec.1")

## perform spatial weighting
prec.1.50 <- SpawExact(precise.data=prec.1,
                      context.id="area.name",
                      contextual.names="prec.1",
                      contextual.weight.matrices=geow,
                      population.weight.names = NULL)

## reshuffle the data and try again
prec.1.shuffle <- prec.1[c(4,2,3,1),]
prec.1.shuffle.50 <- SpawExact(precise.data=prec.1.shuffle,
                               context.id="area.name",
                               contextual.names="prec.1",
                               contextual.weight.matrices=geow,
                               population.weight.names = NULL)
merge.data <- merge(prec.1.50, prec.1.shuffle.50, by="area.name")

## check result correspond to one another
if (!identical(merge.data[["prec.1.1.x"]], merge.data[["prec.1.1.y"]])) {
    stop("order should not be of any importance")
} else {
    print("success.")
}
