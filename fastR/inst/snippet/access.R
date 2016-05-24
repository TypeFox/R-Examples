xm <- matrix(1:16, nrow=4); xm
xm[5]
xm[,2]                   # this is 1 dimensional (a vector)
xm[,2,drop=FALSE]        # this is 2 dimensional (still a matrix)
