require(SwarmSVM)

### Toy Example
local.file.name = tempfile()
download.file("http://www.sfu.ca/~hetongh/data/toydata.RData",local.file.name)
load(local.file.name)
toydata.t = toydata[[2]]
toydata = toydata[[1]]

### ijcnn1

local.file.name = tempfile()
download.file("http://www.sfu.ca/~hetongh/data/ijcnn1.dcsvm.RData",local.file.name)
load(local.file.name)
ijcnn1.t = ijcnn1[[2]]
ijcnn1 = ijcnn1[[1]]


### Covtype
local.file.name = tempfile()
download.file("http://www.sfu.ca/~hetongh/data/covtype.mult.RData",local.file.name)
load(local.file.name)
train.ind = 1:100000
valid.ind = 420001:430000
test.ind = 450001:500000
covtype.train = covtype[train.ind,]
covtype.valid = covtype[valid.ind,]
covtype.test = covtype[test.ind,]

###############
### GaterSVM
###############

### Toy

gaterSVM.model = gaterSVM(x = toydata[,-3], y = toydata[,3], hidden = 10, seed = 0,
                          m = 5, max.iter = 1, learningrate = 0.001, threshold = 0.05,
                          valid.x = toydata.t[,-3], valid.y = toydata.t[,3], verbose = TRUE,
                          stepmax = 100)
table(gaterSVM.model$valid.pred,toydata.t[,3])
#       0    1
# -1 3333    0
# 1     0 6667

### ijcnn1

gaterSVM.model = gaterSVM(x = ijcnn1[,-1], y = ijcnn1[,1], hidden = 50, seed = 0,
                          m = 20, max.iter = 1, learningrate = 0.01, threshold = 0.01,
                          valid.x = ijcnn1.t[,-1], valid.y = ijcnn1.t[,1], verbose = TRUE,
                          stepmax = 100)
table(gaterSVM.model$valid.pred,ijcnn1.t[,1])
#       -1     1
# -1 82081   710
# 1    908  8002
gaterSVM.model$valid.score
# 0.9823556995
gaterSVM.model$time$total.time-gaterSVM.model$time$validation.time
# elapsed 
# 124.547

### covtype

gaterSVM.model = gaterSVM(x = covtype.train[1:50000,-1], y = covtype.train[1:50000,1], 
                          hidden = 50, seed = 0, m = 20, max.iter = 1, 
                          learningrate = 0.03, threshold = 0.05, verbose = TRUE,
                          valid.x = covtype.valid[,-1], 
                          valid.y = covtype.valid[,1],stepmax = 200)
table(gaterSVM.model$valid.pred, covtype.valid[,1])
# pred    0    1
#   -1 5371 2205
#   1   448 1976

gaterSVM.model = gaterSVM(x = covtype.train[,-1], y = covtype.train[,1], 
                          hidden = 150, seed = 0, m = 50, max.iter = 1, 
                          learningrate = 0.03, threshold = 0.05, verbose = TRUE,
                          valid.x = covtype.valid[,-1], 
                          valid.y = covtype.valid[,1],stepmax = 200)
table(gaterSVM.model$valid.pred, covtype.valid[,1])



###############
### LibSVM
###############

LibSVMBenchmark = function(train.x, train.y, test.x, test.y, seed = 1024, ...) {
  set.seed(seed)
  total.time.point = proc.time()
  time.point = proc.time()
  svm.model = svm(x = train.x, y = train.y, ..., probability = FALSE)
  training.time = (proc.time() - time.point)[3]
  preds = as.numeric(predict(svm.model, test.x, probability = FALSE)>0)
  # preds = predict(svm.model, test.x, probability = FALSE)
  print(table(preds, test.y))
  score = sum(diag(table(preds, test.y)))/length(test.y)
  total.time = (proc.time() - total.time.point)[3]
  res = c(score, training.time, total.time)
  names(res) = c("SVM Score", "SVM Training Time", "SVM Total Time")
  return(res)
}

### Toy

LibSVMBenchmark(train.x = toydata[,-3], train.y = toydata[,3], 
                test.x = toydata.t[,-3], test.y = toydata.t[,3])
#       test.y
# preds    0    1
#     0 1536    0
#     1 1797 6667
# SVM Score SVM Training Time    SVM Total Time 
#    0.8203            0.0690            0.2430 

### ijcnn1

LibSVMBenchmark(train.x = ijcnn1[,-1], train.y = ijcnn1[,1], 
                test.x = ijcnn1.t[,-1], test.y = ijcnn1.t[,1],
                gamma = 2, cost = 2)
#      test.y
# preds    -1     1
#     0 82826  2998
#     1   163  5714
#     SVM Score SVM Training Time    SVM Total Time 
#  0.9655292745    862.3340000000    994.8770000000 

### covtype

LibSVMBenchmark(train.x = covtype.train[1:50000,-1], train.y = covtype.train[1:50000,1], 
                test.x = covtype.valid[,-1], test.y = covtype.valid[,1],
                gamma = 2, cost = 2)

LibSVMBenchmark(train.x = covtype.train[,-1], train.y = covtype.train[,1], 
                test.x = covtype.valid[,-1], test.y = covtype.valid[,1],
                gamma = 2, cost = 2)


