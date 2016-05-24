memberships <- function(grid, fact=2){
layers <- layerize(grid)
membership <- aggregate(layers, fact=fact, fun=mean)
return(membership) 
}


