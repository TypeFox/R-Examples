
#f.gaussian.spatial.smoothing.functions.R

#examples of 3 different kernels 
a <- GaussSmoothKernel(voxdim = c(1, 1, 1), ksize = 5, sigma = diag(1, 3))
##f.grid(1,5)
par(mfrow = c(1, 5))
for(i in 1:5){
    image(a[i, , ], zlim = c(min(a), max(a)),axes = FALSE)
    box()
}

a <- GaussSmoothKernel(voxdim = c(1, 1, 1), ksize = 5, sigma = diag(c(0, 2, 2)))
##f.grid(1,5)
par(mfrow = c(1, 5))
for(i in 1:5){
    image(a[i, , ], zlim = c(min(a), max(a)), axes = FALSE)
    box()
}

s <- matrix(2.5, 3, 3)
for(i in 1:3) s[i, i] <- 3
a <- GaussSmoothKernel(voxdim = c(1, 1, 1), ksize = 5, sigma = s)
##f.grid(1,5)
par(mfrow = c(1, 5))
for(i in 1:5){
    image(a[i, , ], zlim = c(min(a), max(a)), axes = FALSE)
    box()
}

#example 1 of GaussSmoothArray

d <- c(10, 10, 10, 20)
mat <- array(rnorm(cumprod(d)[length(d)]), dim = d)
mat[, , 6:10, ] <- mat[, , 6:10 , ] + 3
mask <- array(0, dim = d[1:3])
mask[3:8, 3:8, 3:8] <- 1
b <- GaussSmoothArray(mat, mask = mask, voxdim = c(1, 1, 1), ksize = 5, sigma = diag(1, 3))

f.grid(1,2)
persp(mat[5, , ,1], theta = 90, zlim = c(-2, 6))
persp(b[5, , , 1], theta = 90, zlim = c(-2, 6))         

#example 2 of GaussSmoothArray

d <- c(10, 10, 20, 20)
mat <- array(rnorm(cumprod(d)[length(d)]), dim = d)
mask <- array(1, dim = d[1:3])
#mask[2:9,2:9,2:19]<-1
b <- GaussSmoothArray(mat, voxdim = c(1, 1, 1), ksize = 5, sigma = diag(1, 3), mask = mask, var.norm = TRUE)

##f.grid(1,2)
par(mfrow = c(1, 2))
persp(mat[5, , , 1], theta = 90, zlim = c(-2, 6))
persp(b[5, , , 1], theta = 90, zlim = c(-2, 6))         

#non-linear smoothing

#3D array
d <- rep(10, 3)
a <- array(3, dim = d)
a[, 5:10, 5:10] <- 7
a <- a + array(rnorm(n = 1000, sd = 1), dim = d)

h <- NonLinearSmoothArray(a, voxdim = c(1, 1, 1), radius = 2, sm = 2)

par(mfrow = c(2, 2))
image(a[1, , ], zlim = c(-1, 12)) ; title("Before smoothing")
image(h[1, , ], zlim = c(-1, 12)) ; title("After smoothing")
persp(a[1, , ], zlim = c(-1, 12))
persp(h[1, , ], zlim = c(-1, 12))

#4D array
d <- c(10, 10, 10, 20)
a <- array(1, dim = d)
a[, , 6:10, ] <- 2
a <- a + array(rnorm(20000, sd = .1), dim = d)

h <- NonLinearSmoothArray(a, voxdim = c(1, 1, 1), radius = 2, sm = 3)

f.grid(10,10)
for(i in 1:10){
    for(j in 10:1){
        plot(a[1, i, j, ], type = "l", ylim = c(0, 3), axes = FALSE) ; box()
        lines(h[1, i, j, ], col = 2)
    }
}
