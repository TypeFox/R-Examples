real.dist <-
function(samples, grid, 
                      samples.names = rownames(samples),
                      grid.names = rownames(grid)) {

    d.real <- matrix(0, nrow(grid), nrow(samples))
    colnames(d.real) <- samples.names
    rownames(d.real) <- grid.names

    for (i in 1:nrow(samples)) {
        local <- as.double(samples[i,])
        d.real[,i] <- ((grid[,1]-local[1])**2 + (grid[,2]-local[2])**2)**0.5
    }

    return(d.real)
}
