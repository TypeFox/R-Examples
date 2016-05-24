star3d <-
function (data, vars = NULL, scale = 1) 
{
#    require(rgl)
    p = dim(data)[2] - 1
    n = dim(data)[1]
    theta = 0
    phi = 0
    dtheta = 2 * pi/p
    dphi = pi/p
    u = matrix(0, p, 3)
    for (i in 1:p) {
        u[i, 1] = cos(theta) * sin(phi)
        u[i, 2] = sin(theta) * sin(phi)
        u[i, 3] = cos(phi)
        theta = theta + dtheta
        phi = phi + dphi
    }
    for (j in vars) u[j, ] = u[j, ] * scale
    u1 = rbind(c(0, 0, 0), u)
    #print(u1)
 #   library(rgl)
    rgl::rgl.open()
    u1 = t(u1)
    #print(u1)
    x = u1[1, ]
    y = u1[2, ]
    z = u1[3, ]
    j1 = rep(1, p)
    j2 = 2:(p + 1)
    j12 = cbind(j1, j2)
    j = as.vector(t(j12))
    labels = c("o", paste("V", 1:p, sep = ""))
    maxval = apply(data[, 1:p], 2, max)
    minval = apply(data[, 1:p], 2, min)
    rango = maxval - minval
    mat1 = scale(data[, 1:p], center = minval, scale = rango)
    mat1 = mat1 %*% u
    rgl::plot3d(mat1, col = as.numeric(factor(data[, p + 1])) + 1, add = TRUE, size = 5)
    rgl::text3d(x, y, z, labels)
    rgl::segments3d(x[j], y[j], z[j])
    rgl::identify3d(mat1,labels=1:n)
    rgl::rgl.close()
}
