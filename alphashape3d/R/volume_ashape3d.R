volume_ashape3d <-
function (as3d, byComponents = FALSE, indexAlpha = 1) 
{
    if (class(as3d) != "ashape3d") {
        cat("Argument is not of class ashape3d.\n")
        return(invisible())
    }
    tetra = as3d$tetra
    x <- as3d$x
    if (class(indexAlpha) == "character") 
        if (indexAlpha == "ALL" | indexAlpha == "all") 
            indexAlpha = 1:length(as3d$alpha)
    if (any(indexAlpha > length(as3d$alpha)) | any(indexAlpha <= 
        0)) {
        if (max(indexAlpha) > length(as3d$alpha)) 
            error = max(indexAlpha)
        else error = min(indexAlpha)
        stop(paste("indexAlpha out of bound : valid range = 1:", 
            length(as3d$alpha), ", problematic value = ", error, 
            sep = ""), call. = TRUE)
    }
    for (ii in 1:dim(tetra)[1]) {
        x1 = x[tetra[ii, 1], ]
        x2 = x[tetra[ii, 2], ]
        x3 = x[tetra[ii, 3], ]
        x4 = x[tetra[ii, 4], ]
        tetra[ii, 5] = abs((x4[1] - x1[1]) * ((x2[2] - x1[2]) * 
            (x3[3] - x1[3]) - (x2[3] - x1[3]) * (x3[2] - x1[2])) + 
            (x4[2] - x1[2]) * ((x2[3] - x1[3]) * (x3[1] - x1[1]) - 
                (x2[1] - x1[1]) * (x3[3] - x1[3])) + (x4[3] - 
            x1[3]) * ((x2[1] - x1[1]) * (x3[2] - x1[2]) - (x2[2] - 
            x1[2]) * (x3[1] - x1[1])))/6
    }
    volume = NULL
    if (byComponents) {
        components = components_ashape3d(as3d, indexAlpha = indexAlpha)
        indexComponents = 0
        if (length(indexAlpha) == 1) 
            components = list(components)
        for (iAlpha in indexAlpha) {
            indexComponents = indexComponents + 1
            volumeComponents = rep(0, max(components[[indexComponents]]))
            for (ii in 1:dim(tetra)[1]) {
                if (tetra[ii, 5 + iAlpha] == 1) 
                  volumeComponents[components[[indexComponents]][tetra[ii, 
                    1]]] = volumeComponents[components[[indexComponents]][tetra[ii, 
                    1]]] + tetra[ii, 5]
            }
            if (length(indexAlpha) > 1) {
                volume = c(volume, list(volumeComponents))
            }
            else {
                volume = volumeComponents
            }
        }
    }
    else {
        for (iAlpha in indexAlpha) volume = c(volume, sum(tetra[tetra[, 
            5 + iAlpha] == 1, 5]))
    }
    return(volume)
}
