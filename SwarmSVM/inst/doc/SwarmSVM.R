## ---- message=FALSE------------------------------------------------------
require(SwarmSVM)
data(svmguide1)

## ------------------------------------------------------------------------
head(svmguide1[[1]])

## ------------------------------------------------------------------------
head(svmguide1[[2]])

## ------------------------------------------------------------------------
svmguide1.t = svmguide1[[2]]
svmguide1 = svmguide1[[1]]

## ------------------------------------------------------------------------
csvm.obj = clusterSVM(x = svmguide1[,-1], y = svmguide1[,1], type = 1,
                      valid.x = svmguide1.t[,-1],valid.y = svmguide1.t[,1], 
                      seed = 1, verbose = 1, centers = 8)
csvm.obj$valid.score

## ------------------------------------------------------------------------
class(svmguide1)

## ------------------------------------------------------------------------
csvm.obj = clusterSVM(x = as.matrix(svmguide1[,-1]), y = svmguide1[,1], type = 1,
                      valid.x = as.matrix(svmguide1.t[,-1]),valid.y = svmguide1.t[,1], 
                      seed = 1, verbose = 1, centers = 8)
csvm.obj$valid.score

## ------------------------------------------------------------------------
cluster.fun = stats::kmeans

cluster.predict = function(x, cluster.object) {
  centers = cluster.object$centers
  eucliDist = function(x, centers) apply(centers, 1, function(C) colSums( (t(x)-C)^2 ))
  euclidean.dist = eucliDist(x, centers)
  result = max.col(-euclidean.dist)
  return(result)
}

## ------------------------------------------------------------------------
csvm.obj = clusterSVM(x = svmguide1[,-1], y = svmguide1[,1], centers = 8, seed = 1,
                      cluster.fun = cluster.fun, cluster.predict = cluster.predict,
                      valid.x = svmguide1.t[,-1], valid.y = svmguide1.t[,1])
csvm.obj$valid.score

## ------------------------------------------------------------------------
dcsvm.model = dcSVM(x = svmguide1[,-1], y = svmguide1[,1],
                    k = 4, max.levels = 4, seed = 0, cost = 32, gamma = 2,
                    kernel = 3,early = 0, m = 800,
                    valid.x = svmguide1.t[,-1], valid.y = svmguide1.t[,1])
dcsvm.model$valid.score

## ------------------------------------------------------------------------
dcsvm.model = dcSVM(x = as.matrix(svmguide1[,-1]), y = svmguide1[,1], 
                    k = 10, max.levels = 1, 
                    early = 1, gamma = 2, cost = 32, tolerance = 1e-2, m = 800, 
                    valid.x = svmguide1.t[,-1], valid.y = svmguide1.t[,1])
dcsvm.model$valid.score
dcsvm.model$time$total.time

## ------------------------------------------------------------------------
dcsvm.model = dcSVM(x = as.matrix(svmguide1[,-1]), y = svmguide1[,1], 
                    k = 10, max.levels = 1, 
                    early = 1, gamma = 2, cost = 32, tolerance = 1e-2, m = 800, 
                    valid.x = svmguide1.t[,-1], valid.y = svmguide1.t[,1])
dcsvm.model$valid.score
dcsvm.model$time$total.time

## ------------------------------------------------------------------------
local.file.name = tempfile()
download.file("http://www.sfu.ca/~hetongh/data/ijcnn1.dcsvm.RData",local.file.name)
load(local.file.name)
ijcnn1.t = ijcnn1[[2]]
ijcnn1 = ijcnn1[[1]]

## ------------------------------------------------------------------------
dcsvm.model = dcSVM(x = ijcnn1[,-1], y = ijcnn1[,1], k = 10, max.levels = 1, seed = 1024,
                    early = 1, gamma = 2, cost = 32, tolerance = 1e-2, m = 5000, scale = FALSE,
                    valid.x = ijcnn1.t[,-1], valid.y = ijcnn1.t[,1])
dcsvm.model$valid.score
dcsvm.model$time$total.time

## ------------------------------------------------------------------------
dcsvm.model = dcSVM(x = ijcnn1[,-1], y = ijcnn1[,1], k = 10, max.levels = 1, seed = 1024,
                    early = 1, gamma = 2, cost = 32, tolerance = 1e-2, m = 5000, scale = TRUE,
                    valid.x = ijcnn1.t[,-1], valid.y = ijcnn1.t[,1])
dcsvm.model$valid.score
dcsvm.model$time$total.time

## ------------------------------------------------------------------------
gaterSVM.model = gaterSVM(x = svmguide1[,-1], y = svmguide1[,1], hidden = 10, seed = 0,
                          m = 10, max.iter = 3, learningrate = 0.01, threshold = 1, stepmax = 1000,
                          valid.x = svmguide1.t[,-1], valid.y = svmguide1.t[,1], verbose = TRUE)
gaterSVM.model$valid.score

