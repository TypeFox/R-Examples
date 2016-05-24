cluster.proportion <- function(result, solution, new.order){

    regions <- matrix(NA, nrow=length(result$na.rows), ncol=solution)
    if(result$use.som){
	regions[!result$na.rows,] <- som.regions.for.values(result=result, solution=solution, new.order, get.membership=TRUE)
    } else {
       regions[!result$na.rows,] <- change.order.clusters(result$cluster.assignment[[solution]], new.order)$membership
    }


    if(result$multi.model){
        regions <- drop(add.dimension(regions,t.o=result))
        region.proportion <- apply(regions, MARGIN=2, FUN=colSums,na.rm=TRUE)/NROW(regions)
    } else {
        region.proportion <- t(regions)
    }

    
    return(region.proportion)

}
