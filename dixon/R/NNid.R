`NNid` <-
function (xy, splancs = TRUE) 
{
    if (splancs) {
        nndistG(xy)$neighs
    }
##    else {
##        temp <- #####find.####neighbor(xy, k = 2)
##        as.vector(temp[temp[, 1] != temp[, 2], 2])
##    }
}

