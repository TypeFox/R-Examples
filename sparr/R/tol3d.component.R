
tol3d.component <- function(locations,gridx,gridy,gridz,raise){
    nearestZ <- c()
    for(i in 1:nrow(locations)){
        nearestZ <- append(nearestZ,nnwhich(X=c(locations[i,1],gridx),Y=c(locations[i,2],gridy))[1]-1)
    }
    return(gridz[nearestZ]+raise)
}
