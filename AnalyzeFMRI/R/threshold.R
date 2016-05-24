## Functions for simulating and thresholding Gaussian Random Fields 

Sim.3D.GRF <- function(d, voxdim, sigma, ksize, mask = NULL, type = c("field", "max")) {

    ## simulates a GRF with covariance matrix sigma of dimension d with voxel dimensions voxdim
    ## if type = "max" then just the maximum of the field is returned
    ## if type = "filed" then just the filed AND the maximum of the field are returned
    
    if(!is.null(mask) && sum(d[1:3] == dim(mask)) < 3)  return("mask is wrong size")
    if(length(d) != 3 && length(d) != 4) return("array x should be 3D or 4D")
    if(length(d) == 3) {
        d <- c(d, 1)
        tmp <- 1
    }
    
    if((2 * floor(ksize / 2)) == ksize) stop(paste("ksize must be odd"))
    
    if(is.null(mask)) mask <- array(1, dim = d[1:3])
    space <- 1 + (type == "field") * prod(d)
    
    filtermat <- GaussSmoothKernel(voxdim, ksize, sigma)
    
    a <- .C("sim_grf",
            as.integer(d),
            as.double(aperm(filtermat, c(3, 2, 1))),
            as.integer(ksize),
            as.integer(aperm(mask, c(3, 2, 1))),
            as.integer((type == "field")),
            mat = double(space),
            max = double(1),
            PACKAGE = "AnalyzeFMRI")
    
    if(type == "field") {
        mat <- array(a$mat, dim = d[4:1])
        mat <- aperm(mat, 4:1)
        if(tmp == 1) mat <- mat[, , , 1]
        return(list(mat = mat, max = a$max))
    }
    
    return(list(mat = NULL, max = a$max))
    
    
}


SmoothEst <- function (mat, mask, voxdim, method = "Forman")
{
  ## Estimate the variance-covariance matrix of a Gaussian random field

  ## set-up ## 
  x <- dim(mat)[1]
  y <- dim(mat)[2]
  z <- dim(mat)[3]
  b1 <- array(0, dim = dim(mat) + 2)
  b2 <- array(0, dim = dim(mat) + 2)
  m1 <- array(0, dim = dim(mat) + 2)
  m2 <- array(0, dim = dim(mat) + 2)
  
  ## X ##
  b1[2:(x + 1), 2:(y + 1), 2:(z + 1)] <- mat
  b2[1:(x)    , 2:(y + 1), 2:(z + 1)] <- mat
  m1[2:(x + 1), 2:(y + 1), 2:(z + 1)] <- mask
  m2[1:(x)    , 2:(y + 1), 2:(z + 1)] <- mask
  m3 <- (m1 + m2) == 2
  if(method == "Forman") {
      x.v0 <- mean((b1[m1 == 1])^2)
      x.v1 <- mean(((b1 - b2)[m3 == 1])^2)
      xx <- -(voxdim[1]^2) / (4 * log(1 - x.v1 / (2 * x.v0)))
  } 
  if(method == "Friston") 
      xx <- mean(((b1 - b2)[m3 == 1])^2) / (voxdim[1]^2)
  
  
  ## Y ##
  b1[, , ] <- 0
  b2[, , ] <- 0
  m1[, , ] <- 0
  m2[, , ] <- 0
  b1[2:(x + 1), 2:(y + 1), 2:(z + 1)] <- mat
  b2[2:(x + 1), 1:(y)    , 2:(z + 1)] <- mat
  m1[2:(x + 1), 2:(y + 1), 2:(z + 1)] <- mask
  m2[2:(x + 1), 1:(y)    , 2:(z + 1)] <- mask
  m3 <- (m1 + m2) == 2
  if(method == "Forman") {
    y.v0 <- mean((b1[m1 == 1])^2)
    y.v1 <- mean(((b1 - b2)[m3 == 1])^2)
    yy <- -(voxdim[2]^2) / (4 * log(1 - y.v1 / (2 * y.v0)))
  }
  if(method == "Friston")
      yy <- mean(((b1 - b2)[m3 == 1])^2) / (voxdim[2]^2)
  
  ## Z ##
  b1[, , ] <- 0
  b2[, , ] <- 0
  m1[, , ] <- 0
  m2[, , ] <- 0
  b1[2:(x + 1), 2:(y + 1), 2:(z + 1)] <- mat
  b2[2:(x + 1), 2:(y + 1), 1:(z)    ] <- mat
  m1[2:(x + 1), 2:(y + 1), 2:(z + 1)] <- mask
  m2[2:(x + 1), 2:(y + 1), 1:(z)    ] <- mask
  m3 <- (m1 + m2) == 2
  if(method == "Forman") {
      z.v0 <- mean((b1[m1 == 1])^2)
      z.v1 <- mean(((b1 - b2)[m3 == 1])^2)
      zz <- -(voxdim[3]^2) / (4 * log(1 - z.v1 / (2 * z.v0)))
    }
  if(method == "Friston") 
      zz <- mean(((b1 - b2)[m3 == 1])^2) / (voxdim[3]^2)
  
  ## output ##
  if(method == "Forman") 
    sigma <- diag(c(xx, yy, zz))
  if(method == "Friston")
    sigma <- solve(diag(2 * c(xx, yy, zz)))
  
  return(sigma)
  
}


Threshold.Bonferroni <- function(p.val, n, type = c("Normal", "t", "F"), df1 = NULL, df2 = NULL) {

    ## calculate the Bonferroni threshold for n iid tests to give a p-value of p.val
    ## type specifies the univariate distribution of the test statistics under consideration
    
    if(type == "Normal") return(qnorm(1 - p.val / n))
    if(type == "t") return(qt(1 - p.val / n, df = df1))
    if(type == "F") return(qf(1 - p.val / n, df1 = df1, df2 = df2))

}

EC.3D <- function(u, sigma, voxdim = c(1, 1, 1), num.vox, type = c("Normal", "t"), df = NULL) {

    ## The Expectation of the Euler Characteristic for a 3D Random Field above a threshold u
    ## type specifies the marginal distribution of the field
    
    V <- prod(voxdim) * num.vox
    S <- sqrt(det(solve(2 * sigma)))
    EC <- switch(type[1],
                 Normal = V * S * (u^2 - 1) * exp(-u^2 / 2) / (4 * pi^2),
                 t = V * S * (1 + (u^2 / df))^(-(df - 1) / 2) * (u^2 * (df - 1) / df - 1) / (4 * pi^2)
                 )
    return(EC)
}


Threshold.RF <- function(p.val, sigma, voxdim = c(1, 1, 1), num.vox, type = c("Normal", "t"), df = NULL) {

    ## calculates the Random Field theory threshold to give a p-value of p.val
    ## type specifies the marginal distribution of the field
    
    EC.func <- function(u , sigma, voxdim, num.vox, p.val) (EC.3D(u, sigma, voxdim, num.vox, type[1], df) - p.val)^2

    threshold <- optimize(f = EC.func, interval = c(2,10), sigma = sigma, voxdim = voxdim, num.vox = num.vox, p.val = p.val, tol = 1e-9)$minimum
    
    return(threshold)
}


Threshold.FDR <- function(x, q, cV.type = 2, type = c("Normal", "t", "F"), df1 = NULL, df2 = NULL) {

    ## calculates the FDR threshold for a vector of p-values in x
    ## q specifies the desired FDR
    ## cV specfies the type of FDR threshold used (See Genovese et al. (2002))
    
    if(type == "Normal") p <- sort(1 - pnorm(x))
    if(type == "t") p <- sort(1 - pt(x, df = df1))
    if(type == "F") p <- sort(1 - pf(x, df1 = df1, df2 = df2))

    V <- length(p)
    cV <- switch(cV.type, 1, log(V) + 0.5772)
    
    i <- 1
    while (p[i] <= (i * q) / (V * cV)) i <- i + 1
    i <- max(i - 1, 1)

    if(type == "Normal") thr <- qnorm(1 - p[i])
    if(type == "t") thr <- qt(1 - p[i], df = df1)
    if(type == "F") thr <- qf(1 - p[i], df1 = df1, df2 = df2)

    return(thr)

}

Sim.3D.GammaRF <- function(d, voxdim, sigma, ksize, mask, shape, rate) {

  ## simulates a smooth Gamma distributed random field by simulating a GRF and
  ## transforming each statistic value to be a Gammma
  
  field <- Sim.3D.GRF(d = d, voxdim = voxdim, sigma = sigma, ksize = 9,
                      mask = mask, type = "field")$mat

  gamma.field <- qgamma(pnorm(field), shape = shape, rate = rate)
  gamma.field <- gamma.field * mask

  return(gamma.field)

}
