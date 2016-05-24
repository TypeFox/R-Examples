## -----------------------------------
## Simulate Tumor Data and write a csv file
## -----------------------------------

simtum = function(num.tum = 3, width = 70, height = 70, radius.l = 4, radius.u = 12, seed.center = NULL){
    radius = runif(num.tum,radius.l, radius.u) # Radius of the tumors
    tum = data.frame(row = rep(1:height, each = width), col = rep(1:height, width)) 
    tum$pixel = rnorm(width * height, 0, 0.5) # White noises
    tum$is = FALSE # Is it a tumor cell?
    if(!missing(seed.center)) set.seed(seed.center)
    centers = list(row = c(sample(1:height, num.tum, replace = TRUE)),
        col = c(sample(1:width, num.tum, replace = TRUE))) # Centers of the tumors
                                        # Generate pixels for tumors.
    for(i in 1:num.tum){
        in.circle = which(sqrt((tum$row-centers$row[i])^2 + (tum$col-centers$col[i])^2) <= radius[i])
        tum$pixel[in.circle] = tum$pixel[in.circle] + rnorm(length(in.circle),4,1)
        tum$is[in.circle] = TRUE # Is the pixel a tumor pixel?
    } # Generate tumor pixels
    tum.range = range(tum$pixel)
    tum$pixel = tum$pixel - tum.range[1]
    tum$pixel = tum$pixel / (tum.range[2] - tum.range[1])
    return(invisible(tum$pixel))
}

## image(matrix(simtum(),70))
    
## ## -----------------------------------
## ## Simulate Tumor with different tumor locations and outlier locations.
## ## -----------------------------------
## tum.multiple = replicate(1000, simtum())
## tum.multiple.o = tum.multiple

## add.outliers = function(x, value = 1, per = 0.05){
##     x[sample(1:length(x), (per * length(x)))] = value
##     return(x)
## }

## tum.multiple = apply(tum.multiple.o,2,add.outliers)

## tumors = list(NA, 1000)
## for(i in 1:1000){
##     tumors[[i]] = list(original = matrix(tum.multiple.o[,i],70,70), corrupt = matrix(tum.multiple[,i],70,70))
## }

## save.image("tumors.Rdata")




