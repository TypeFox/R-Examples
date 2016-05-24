dist.to.knn <-
function (dataset, neighbors) 
{
    numrow = dim(dataset)[1]
    knndist = rep(0, 0)
    for (i in 1:numrow) {
        neighdist = knneigh.vect(dataset[i, ], dataset, neighbors)
        if (i == 2) {
            if (length(knndist) < length(neighdist)) {
                z = length(neighdist) - length(knndist)
                zeros = rep(0, z)
                knndist = c(knndist, zeros)
            }
            else if (length(knndist) > length(neighdist)) {
                z = length(knndist) - length(neighdist)
                zeros = rep(0, z)
                neighdist = c(neighdist, zeros)
            }
        }
        else {
            if (i != 1) {
                if (dim(knndist)[1] < length(neighdist)) {
                  z = (length(neighdist) - dim(knndist)[1])
                  zeros = rep(0, z * dim(knndist)[2])
                  zeros = matrix(zeros, z, dim(knndist)[2])
                  knndist = rbind(knndist, zeros)
                }
                else if (dim(knndist)[1] > length(neighdist)) {
                  z = (dim(knndist)[1] - length(neighdist))
                  zeros = rep(0, z)
                  neighdist = c(neighdist, zeros)
                }
            }
        }
        knndist = cbind(knndist, neighdist)
    }
    return(knndist)
}
