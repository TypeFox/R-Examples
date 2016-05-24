require(SwarmSVM)

## SVMGUIDE1
local.file.name = tempfile()
download.file("http://www.sfu.ca/~hetongh/data/svmguide1.RData",local.file.name)
load(local.file.name)

svmguide1.t = svmguide1[[2]]
svmguide1 = svmguide1[[1]]
csvm.obj = clusterSVM(x = svmguide1[,-1], y = svmguide1[,1],
                      centers = 8, iter.max = 1000, seed = 512, verbose = 0,
                      valid.x = svmguide1.t[,-1],valid.y = svmguide1.t[,1])
# Time for Clustering: 0.022 secs
# Time for Transforming: 0.026 secs
# Time for Liblinear: 0.202 secs
# Time for Validation: 0.049 secs
# 
# Total Time: 0.299 secs
# Accuracy Score: 0.80425 

# Dense Matrix
csvm.obj = clusterSVM(x = svmguide1[,-1], y = svmguide1[,1], sparse = FALSE,
                      centers = 8, iter.max = 1000, seed = 512,
                      valid.x = svmguide1.t[,-1],valid.y = svmguide1.t[,1])
# Time for Clustering: 0.022 secs
# Time for Transforming: 0.003 secs
# Time for Liblinear: 0.128 secs
# Time for Validation: 0.018 secs
# 
# Total Time: 0.173 secs
# Accuracy Score: 0.80425

# Self-Defined clustering algorithm
cluster.fun = function(x, centers, ...) {
  x = as.matrix(x)
  kernl.result = kernlab::kkmeans(x, centers, ...)
  result = list()
  result$cluster = kernl.result@.Data
  result$centers = kernl.result@centers
  return(result)
}
csvm.obj = clusterSVM(x = svmguide1[,-1], y = svmguide1[,1], seed = 512,
                      cluster.fun = cluster.fun, centers = 8, 
                      valid.x = svmguide1.t[,-1],valid.y = svmguide1.t[,1])
# Time for Clustering: 11.925 secs
# Time for Transforming: 0.023 secs
# Time for Liblinear: 0.144 secs
# Time for Validation: 0.034 secs
# 
# Total Time: 12.127 secs
# Accuracy Score: 0.81325 

# Using stats::kmeans
csvm.obj = clusterSVM(x = svmguide1[,-1], y = svmguide1[,1], seed = 512,
                      cluster.fun = stats::kmeans, centers = 8, 
                      valid.x = svmguide1.t[,-1],valid.y = svmguide1.t[,1])
# KMeans::Cluster(): converged after 47 iterations.
# exiting from: cluster.fun(x, ...)
# Time for Clustering: 20.174 secs
# Time for Transforming: 12.39 secs
# Time for Liblinear: 12.993 secs
# Time for Validation: 4.106 secs
# 
# Total Time: 49.681 secs
# Accuracy Score: 0.9614 

## IJCNN1
# data(ijcnn1)
local.file.name = tempfile()
download.file("http://www.sfu.ca/~hetongh/data/ijcnn1.RData",local.file.name)
load(local.file.name)

ijcnn1.t = ijcnn1[[2]]
ijcnn1 = ijcnn1[[1]]
csvm.obj = clusterSVM(x = ijcnn1[,-1], y = ijcnn1[,1],
                      centers = 8, iter.max = 1000, seed = 512,
                      valid.x = ijcnn1.t[,-1],valid.y = ijcnn1.t[,1])
# Time for Clustering: 0.2 secs
# Time for Transforming: 0.306 secs
# Time for Liblinear: 1.972 secs
# Time for Validation: 2.059 secs
# 
# Total Time: 4.547 secs
# Accuracy Score: 0.9425742 



## USPS
# data(usps)
local.file.name = tempfile()
download.file("http://www.sfu.ca/~hetongh/data/usps.RData",local.file.name)
load(local.file.name)

usps.t = usps[[2]]
usps = usps[[1]]
csvm.obj = clusterSVM(x = usps[,-1], y = usps[,1],
                      centers = 8, iter.max = 1000, seed = 512,
                      valid.x = usps.t[,-1],valid.y = usps.t[,1])
# Time for Clustering: 1.699 secs
# Time for Transforming: 1.03 secs
# Time for Liblinear: 1.334 secs
# Time for Validation: 0.608 secs
# 
# Total Time: 4.676 secs
# Accuracy Score: 0.9531639



## MNIST
# data(mnist)
local.file.name = tempfile()
download.file("http://www.sfu.ca/~hetongh/data/mnist.RData",local.file.name)
load(local.file.name)

mnist38 = mnist[[1]]
mnist38.t = mnist[[2]]
mnist49 = mnist[[3]]
mnist49.t = mnist[[4]]
mnistoe = mnist[[5]]
mnistoe.t = mnist[[6]]

# 3 vs 8
csvm.obj = clusterSVM(x = mnist38[,-1], y = mnist38[,1],
                      centers = 8, iter.max = 1000, seed = 512,
                      valid.x = mnist38.t[,-1],valid.y = mnist38.t[,1])
# Time for Clustering: 11.777 secs
# Time for Transforming: 2.925 secs
# Time for Liblinear: 2.531 secs
# Time for Validation: 0.74 secs
# 
# Total Time: 17.976 secs
# Accuracy Score: 0.984375 

# 4 vs 9
csvm.obj = clusterSVM(x = mnist49[,-1], y = mnist49[,1],
                      centers = 8, iter.max = 1000, seed = 512,
                      valid.x = mnist49.t[,-1],valid.y = mnist49.t[,1])
# Time for Clustering: 8.846 secs
# Time for Transforming: 2.156 secs
# Time for Liblinear: 2.305 secs
# Time for Validation: 0.73 secs
# 
# Total Time: 14.041 secs
# Accuracy Score: 0.9804119

# O vs E
csvm.obj = clusterSVM(x = mnistoe[,-1], y = mnistoe[,1],
                      centers = 8, iter.max = 1000, seed = 512,
                      valid.x = mnistoe.t[,-1],valid.y = mnistoe.t[,1])
# Time for Clustering: 62.051 secs
# Time for Transforming: 11.932 secs
# Time for Liblinear: 13.056 secs
# Time for Validation: 4.122 secs
# 
# Total Time: 91.175 secs
# Accuracy Score: 0.9611 

