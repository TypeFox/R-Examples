#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%
#%%%%%%%  Assess the validity of the Growing region algorithm
#%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

require(MRIaggr)
require(Rcpp)
optionsMRIaggr(mar = c(3,0,3,0), axes = FALSE, main.legend = "", legend = FALSE, 
               bg = "white", num.main = FALSE, outline.index = TRUE)

# source("E:/Creation_package/Package_MRIaggr/MRIaggr/R/multipar.R")
# source("E:/Creation_package/Package_MRIaggr/MRIaggr/R/Associated_functions.R")
# sourceCpp("E:/Creation_package/Package_MRIaggr/MRIaggr/src/Associated_functions.cpp")

#### 1- simulation study ####


#### settings ####
G <- 2 # nb of groups
n_px <- 30 # width of the area
n <- (n_px * G)^2 # total number of pixel

radius.circle1 <- 15
longueur.circle1 <- 25
largeur.circle1 <- 10

# creation of the field
M_coords <- matrix(1, nrow = n_px * G, ncol = n_px * G)
# coordinates of the pixels
coords <- data.frame(cbind(which(M_coords > 0, arr.ind = TRUE), 1))
names(coords) <- c("i","j","k")
contrast <- rep(0, n)
# creation of the groups
center <- apply(coords, 2, function(x){median(unique(x))})
index.circle1 <- which(sqrt(rowSums(sweep(coords, 2, FUN = "-", center)^2)) < radius.circle1)
index.rectangle1 <-  intersect(which(abs(coords[,1] - center[1]) < longueur.circle1),
                               which(abs(coords[,2] - center[2]) < largeur.circle1))
index.seed <- which(sqrt(rowSums(sweep(coords, 2, FUN="-", center)^2)) < 2)

data.image <- coords
data.image$seed <- 0
data.image$circle1 <- 0
data.image$rectangle1 <- 0
data.image$seed[index.seed] <- 1
data.image$circle1[index.circle1] <- 1
data.image$rectangle1[index.rectangle1] <- 1

data.image$Icircle1_sd01 <- rnorm(n, mean = data.image$circle1, sd = 0.1)
data.image$Irectangle1_sd01 <- rnorm(n, mean = data.image$rectangle1, sd = 0.1)
data.image$IGcircle1_sd01 <- rnorm(n, mean = data.image$circle1, sd = 0.1) - log(sqrt(rowSums(sweep(coords, 2, FUN = "-", center)^2)))
data.image$IGrectangle1_sd01 <- rnorm(n, mean = data.image$rectangle1, sd = 0.1) - log(sqrt(rowSums(sweep(coords, 2, FUN = "-", center)^2)))
data.image$Icircle1_sd025 <- rnorm(n, mean = data.image$circle1, sd = 0.25)
data.image$Irectangle1_sd025 <- rnorm(n, mean = data.image$rectangle1,sd = 0.25)

# display
multiplot(coords, data.image[,"Icircle1_sd01", drop = FALSE])
multiplot(coords, data.image[,"Irectangle1_sd01", drop = FALSE])

multiplot(coords, data.image[,"IGcircle1_sd01", drop = FALSE])
multiplot(coords, data.image[,"IGrectangle1_sd01", drop = FALSE])

multiplot(coords, data.image[,"Icircle1_sd025", drop = FALSE])
multiplot(coords, data.image[,"Irectangle1_sd025", drop = FALSE])

#### neighborhood matrix
W.image <- calcW(coords, 1.875, method = "euclidean", upper = NULL)$W



#### 2- little noise ####

sigma_test <- c(seq(0.05, 0.15, 0.01),seq(0.175, 0.5, 0.025))

#### circle
resSigmaGR1 <- calcSigmaGR(contrast = data.image$Icircle1_sd01, W = W.image, 
                           seed = index.seed, sigma = sigma_test,                      
                           criterion.entropy = TRUE,
                           criterion.Kalinsky = TRUE,
                           criterion.Laboure = TRUE,
                           criterion.transition = TRUE,
                           criterion.sdfront = TRUE)

plotSigmaGR(resSigmaGR1)

multiplot(coords, data.image[,"Icircle1_sd01", drop = FALSE],
          index1 = list(coords = coords[resSigmaGR1$list.GR$Kalinsky,]))
multiplot(coords, data.image[,"Icircle1_sd01", drop = FALSE],
          index1 = list(coords = coords[resSigmaGR1$list.GR$Kalinsky,], outline = FALSE))
multiplot(coords, data.image[,"Icircle1_sd01", drop = FALSE],
          index1 = list(coords = coords[resSigmaGR1$list.GR$entropy,]))


resGR1 <- calcGR(contrast = data.image$Icircle1_sd01, W = W.image, 
                 seed = index.seed, sigma_max = 0.4, iter_max = 40, verbose = TRUE)

multiplot(coords, data.image[,"Icircle1_sd01",drop = FALSE],
          index1 = list(coords = coords[resGR1$GR,], outline = FALSE))

#### rectangle
resSigmaGR1bis <- calcSigmaGR(contrast = data.image$Irectangle1_sd01, W = W.image, 
                           seed = index.seed, sigma = sigma_test,                      
                           criterion.entropy = TRUE,
                           criterion.Kalinsky = TRUE,
                           criterion.Laboure = TRUE,
                           criterion.transition = TRUE,
                           criterion.sdfront = TRUE)

plotSigmaGR(resSigmaGR1bis)

multiplot(coords, data.image[,"Irectangle1_sd01", drop = FALSE],
          index1 = list(coords = coords[resSigmaGR1bis$list.GR$Kalinsky,], outline = TRUE))
multiplot(coords, data.image[,"Irectangle1_sd01", drop = FALSE],
          index1 = list(coords = coords[resSigmaGR1bis$list.GR$entropy,]))

#### 2- moderate noise ####
sigma_test <- c(seq(0.05, 0.15, 0.01),seq(0.175, 0.5, 0.025))

#### circle
resSigmaGR2 <- calcSigmaGR(contrast = data.image$Icircle1_sd025, W = W.image, 
                           seed = index.seed, sigma = sigma_test,                      
                           criterion.entropy = TRUE,
                           criterion.Kalinsky = TRUE,
                           criterion.Laboure = TRUE,
                           criterion.transition = TRUE,
                           criterion.sdfront = TRUE)

plotSigmaGR(resSigmaGR2)

multiplot(coords, data.image[,"Icircle1_sd025", drop = FALSE],
          index1 = list(coords = coords[resSigmaGR2$list.GR$Kalinsky,]))
multiplot(coords, data.image[,"Icircle1_sd025", drop = FALSE],
          index1 = list(coords = coords[resSigmaGR2$list.GR$Kalinsky,], outline = TRUE))
multiplot(coords, data.image[,"Icircle1_sd025", drop = FALSE],
          index1 = list(coords = coords[resSigmaGR2$list.GR$entropy,]))

#### rectangle
resSigmaGR2bis <- calcSigmaGR(contrast = data.image$Irectangle1_sd025, W = W.image, 
                              seed = index.seed, sigma = sigma_test,                      
                              criterion.entropy = TRUE,
                              criterion.Kalinsky = TRUE,
                              criterion.Laboure = TRUE,
                              criterion.transition = TRUE,
                              criterion.sdfront = TRUE)

plotSigmaGR(resSigmaGR2bis)

multiplot(coords, data.image[,"Irectangle1_sd025", drop = FALSE],
          index1 = list(coords = coords[resSigmaGR2bis$list.GR$Kalinsky,], outline = TRUE))
multiplot(coords, data.image[,"Irectangle1_sd025", drop = FALSE],
          index1 = list(coords = coords[resSigmaGR2bis$list.GR$sdfront,], outline = TRUE))
multiplot(coords, data.image[,"Irectangle1_sd025", drop = FALSE],
          index1 = list(coords = coords[resSigmaGR2bis$list.GR$entropy,]))


#### 3- gradiant little noise ####
sigma_test <- c(seq(0.05, 0.15, 0.01),seq(0.175, 0.375, 0.1),seq(0.4, 0.8, 0.025))

#### circle
resSigmaGR3 <- calcSigmaGR(contrast = data.image$IGcircle1_sd01, W = W.image, 
                           seed = index.seed, sigma = sigma_test,                      
                           criterion.entropy = TRUE,
                           criterion.Kalinsky = TRUE,
                           criterion.Laboure = TRUE,
                           criterion.transition = TRUE,
                           criterion.sdfront = TRUE)

plotSigmaGR(resSigmaGR3)

multiplot(coords, data.image[,"IGcircle1_sd01", drop = FALSE],
          index1 = list(coords = coords[resSigmaGR3$list.GR$Kalinsky,]))
multiplot(coords, data.image[,"IGcircle1_sd01", drop = FALSE],
          index1 = list(coords = coords[resSigmaGR3$list.GR$Kalinsky,], outline = TRUE))
multiplot(coords, data.image[,"IGcircle1_sd01", drop = FALSE],
          index1 = list(coords = coords[resSigmaGR3$list.GR$entropy,]))

#### rectangle
resSigmaGR3bis <- calcSigmaGR(contrast = data.image$IGrectangle1_sd01, W = W.image, 
                              seed = index.seed, sigma = sigma_test,                      
                              criterion.entropy = TRUE,
                              criterion.Kalinsky = TRUE,
                              criterion.Laboure = TRUE,
                              criterion.transition = TRUE,
                              criterion.sdfront = TRUE)

plotSigmaGR(resSigmaGR3bis)

multiplot(coords, data.image[,"IGrectangle1_sd01", drop = FALSE],
          index1 = list(coords = coords[resSigmaGR3bis$list.GR$Kalinsky,], outline = TRUE))
multiplot(coords, data.image[,"IGrectangle1_sd01", drop = FALSE],
          index1 = list(coords = coords[resSigmaGR3bis$list.GR$sdfront,], outline = TRUE))
multiplot(coords, data.image[,"IGrectangle1_sd01", drop = FALSE],
          index1 = list(coords = coords[resSigmaGR3bis$list.GR$entropy,], outline = TRUE))

resSigmaGR3bis <- calcSigmaGR(contrast = data.image$IGrectangle1_sd01, W = W.image, 
                              seed = index.seed, sigma = sigma_test, keep.upper = TRUE,                      
                              criterion.entropy = TRUE,
                              criterion.Kalinsky = TRUE,
                              criterion.Laboure = TRUE,
                              criterion.transition = TRUE,
                              criterion.sdfront = TRUE)

multiplot(coords, data.image[,"IGrectangle1_sd01", drop = FALSE],
          index1 = list(coords = coords[resSigmaGR3bis$list.GR$Kalinsky,], outline = FALSE))
multiplot(coords, data.image[,"IGrectangle1_sd01", drop = FALSE],
          index1 = list(coords = coords[resSigmaGR3bis$list.GR$sdfront,], outline = FALSE))
multiplot(coords, data.image[,"IGrectangle1_sd01", drop = FALSE],
          index1 = list(coords = coords[resSigmaGR3bis$list.GR$entropy,], outline = TRUE))

#### 4- initialisation ####