som.regions.for.values <- function(result, solution, new.order=1:solution, get.membership=FALSE){

        clustering <- change.order.clusters(result$cluster.assignment[[solution]], new.order)
        which.som.cell <- with(result$som$visual, 1+x+y*result$som.dim[1])

        if(get.membership){
            clusters <- clustering$membership
            regions <- clusters[which.som.cell,]
        } else {
            clusters <- clustering$cluster
            regions <- clusters[which.som.cell]
        }
        return(regions)
}
