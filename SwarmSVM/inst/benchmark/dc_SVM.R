require(SwarmSVM)
require(e1071)

### ijcnn1

local.file.name = tempfile()
download.file("http://www.sfu.ca/~hetongh/data/ijcnn1.dcsvm.RData",local.file.name)
load(local.file.name)
ijcnn1.t = ijcnn1[[2]]
ijcnn1 = ijcnn1[[1]]

### covtype

local.file.name = tempfile()
download.file("http://www.sfu.ca/~hetongh/data/covtype.RData",local.file.name)
load(local.file.name)
covtype.t = covtype[[2]]
covtype = covtype[[1]]

### Benchmark function

dcSVMBenchmark = function(...) {
  dcsvm.model = dcSVM(..., pre.scale = FALSE, k = 10, max.levels = 1, early = 1)
  early.score = dcsvm.model$valid.score
  early.training = dcsvm.model$time$total.time - dcsvm.model$time$validation.time
  early.time = dcsvm.model$time$total.time
  cat('Early Training Finished.\n')
  dcsvm.model = dcSVM(..., k = 4, pre.scale = FALSE, max.levels = 4, early = 0)
  exact.score = dcsvm.model$valid.score
  exact.training = dcsvm.model$time$total.time - dcsvm.model$time$validation.time
  exact.time = dcsvm.model$time$total.time
  cat('Exact Training Finished.\n')
  res = c(early.score, early.training, early.time, 
          exact.score, exact.training, exact.time)
  names(res) = c("Early Score", "Early Training", "Early Time", 
                 "Exact Score","Exact Training", "Exact Time")
  return(res)
}

##############
#### DC SVM
##############

### ijcnn1

dcSVMBenchmark(x = ijcnn1[,-1], y = ijcnn1[,1], seed = 1024,
               gamma = 2^6, cost = 2^(-10), tolerance = 1e-2, m = 2000, scale = FALSE,
               valid.x = ijcnn1.t[,-1], valid.y = ijcnn1.t[,1])
# Early Score Early Training     Early Time    Exact Score Exact Training     Exact Time 
# 0.9370127     79.8400000    122.7020000      0.9371544   1202.9220000   1582.5420000

dcSVMBenchmark(x = ijcnn1[,-1], y = ijcnn1[,1], seed = 1024,
               gamma = 2^10, cost = 2^(-10), tolerance = 1e-2, m = 2000, scale = FALSE,
               valid.x = ijcnn1.t[,-1], valid.y = ijcnn1.t[,1])
# Early Score Early Training     Early Time    Exact Score Exact Training     Exact Time 
# 0.9703929     87.9920000    131.7000000      0.9704801    983.8050000   1344.1000000

dcSVMBenchmark(x = ijcnn1[,-1], y = ijcnn1[,1], seed = 1024,
               gamma = 2^10, cost = 2^(-6), tolerance = 1e-2, m = 2000, scale = FALSE,
               valid.x = ijcnn1.t[,-1], valid.y = ijcnn1.t[,1])
# Early Score Early Training     Early Time    Exact Score Exact Training     Exact Time 
# 0.9704474     89.2600000    133.4160000      0.9705347   1022.8440000   1384.3850000 

dcSVMBenchmark(x = ijcnn1[,-1], y = ijcnn1[,1], seed = 1024,
               gamma = 2, cost = 2, tolerance = 1e-2, m = 2000, scale = FALSE,
               valid.x = ijcnn1.t[,-1], valid.y = ijcnn1.t[,1])
# Early Score Early Training     Early Time    Exact Score Exact Training     Exact Time 
# 0.9879500      8.6720000     15.2110000      0.9883753     41.9140000     78.0440000 

dcSVMBenchmark(x = ijcnn1[,-1], y = ijcnn1[,1], seed = 1024,
               gamma = 2, cost = 2^6, tolerance = 1e-2, m = 5000, scale = FALSE,
               valid.x = ijcnn1.t[,-1], valid.y = ijcnn1.t[,1])
# Early Score Early Training     Early Time    Exact Score Exact Training     Exact Time 
# 0.9834462      7.0340000     11.4670000      0.9833044     25.4920000     44.6110000  

dcSVMBenchmark(x = ijcnn1[,-1], y = ijcnn1[,1], seed = 1024,
               gamma = 2^6, cost = 2^10, tolerance = 1e-2, m = 5000, scale = FALSE,
               valid.x = ijcnn1.t[,-1], valid.y = ijcnn1.t[,1])
# Early Score Early Training     Early Time    Exact Score Exact Training     Exact Time 
# 0.9816142     71.1820000     97.0980000      0.9817232    446.3780000    668.8540000

### covtype

dcSVMBenchmark(x = covtype[,-1], y = covtype[,1], seed = 1024,
               gamma = 2^10, cost = 2^(-10), tolerance = 1e-2, m = 2000, scale = FALSE,
               valid.x = covtype.t[,-1], valid.y = covtype.t[,1])

dcSVMBenchmark(x = covtype[,-1], y = covtype[,1], seed = 1024,
               gamma = 2^10, cost = 2^(-6), tolerance = 1e-2, m = 2000, scale = FALSE,
               valid.x = covtype.t[,-1], valid.y = covtype.t[,1])

dcSVMBenchmark(x = covtype[,-1], y = covtype[,1], seed = 1024,
               gamma = 2^6, cost = 2, tolerance = 1e-2, m = 2000, scale = FALSE,
               valid.x = covtype.t[,-1], valid.y = covtype.t[,1])


dcSVMBenchmark(x = covtype[,-1], y = covtype[,1], seed = 1024,
               gamma = 2^6, cost = 2^6, tolerance = 1e-2, m = 5000, scale = FALSE,
               valid.x = covtype.t[,-1], valid.y = covtype.t[,1])


dcSVMBenchmark(x = covtype[,-1], y = covtype[,1], seed = 1024,
               gamma = 2^6, cost = 2^10, tolerance = 1e-2, m = 5000, scale = FALSE,
               valid.x = covtype.t[,-1], valid.y = covtype.t[,1])

################
#### SVM
################

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


### ijcnn1

LibSVMBenchmark(train.x = ijcnn1[,-1], train.y = ijcnn1[,1], 
                test.x = ijcnn1.t[,-1], test.y = ijcnn1.t[,1],
                gamma = 2^6, cost = 2^(-10), scale = FALSE)
# SVM Score SVM Training Time    SVM Total Time 
# 0.9049956       762.3280000       843.3250000 

LibSVMBenchmark(train.x = ijcnn1[,-1], train.y = ijcnn1[,1], 
                test.x = ijcnn1.t[,-1], test.y = ijcnn1.t[,1],
                gamma = 2^10, cost = 2^(-10))
# SVM Score SVM Training Time    SVM Total Time 
# 0.9049956       125.3940000       191.7660000 

LibSVMBenchmark(train.x = ijcnn1[,-1], train.y = ijcnn1[,1], 
                test.x = ijcnn1.t[,-1], test.y = ijcnn1.t[,1],
                gamma = 2^10, cost = 2^(-6))
# SVM Score SVM Training Time    SVM Total Time 
# 0.9049956      1122.2480000      1433.9050000 

LibSVMBenchmark(train.x = ijcnn1[,-1], train.y = ijcnn1[,1], 
                test.x = ijcnn1.t[,-1], test.y = ijcnn1.t[,1],
                gamma = 2, cost = 2)
# SVM Score SVM Training Time    SVM Total Time 
# 0.9655293       844.0770000       969.9170000 

LibSVMBenchmark(train.x = ijcnn1[,-1], train.y = ijcnn1[,1], 
                test.x = ijcnn1.t[,-1], test.y = ijcnn1.t[,1],
                gamma = 2, cost = 2^6)
# SVM Score SVM Training Time    SVM Total Time 
# 9.681683e-01      1.507333e+04      1.520635e+04 

LibSVMBenchmark(train.x = ijcnn1[,-1], train.y = ijcnn1[,1], 
                test.x = ijcnn1.t[,-1], test.y = ijcnn1.t[,1],
                gamma = 2^6, cost = 2^10)
# SVM Score SVM Training Time    SVM Total Time 
# 0.9565108      2417.9580000      2640.9240000 

### covtype

LibSVMBenchmark(train.x = covtype[,-1], train.y = covtype[,1], 
                test.x = covtype.t[,-1], test.y = covtype.t[,1],
                gamma = 2^6, cost = 2^(-10))

LibSVMBenchmark(train.x = covtype[,-1], train.y = covtype[,1], 
                test.x = covtype.t[,-1], test.y = covtype.t[,1],
                gamma = 2^10, cost = 2^(-10))
# Early Score Early Training     Early Time    Exact Score Exact Training     Exact Time 
# 9.554569e-01   7.602972e+03   8.101556e+03   9.554741e-01   7.071668e+04   7.435093e+04

LibSVMBenchmark(train.x = covtype[,-1], train.y = covtype[,1], 
                test.x = covtype.t[,-1], test.y = covtype.t[,1],
                gamma = 2^10, cost = 2^(-6))

LibSVMBenchmark(train.x = covtype[,-1], train.y = covtype[,1], 
                test.x = covtype.t[,-1], test.y = covtype.t[,1],
                gamma = 2, cost = 2)

LibSVMBenchmark(train.x = covtype[,-1], train.y = covtype[,1], 
                test.x = covtype.t[,-1], test.y = covtype.t[,1],
                gamma = 2, cost = 2^6)
# Early Score Early Training     Early Time    Exact Score Exact Training     Exact Time 
# 0.9834462      7.0340000     11.4670000      0.9833044     25.4920000     44.6110000 

LibSVMBenchmark(train.x = covtype[,-1], train.y = covtype[,1], 
                test.x = covtype.t[,-1], test.y = covtype.t[,1],
                gamma = 2^6, cost = 2^10)
