#*# This demo shows a simple data mining process (phase 1 of TDMR) for classification on
#*# dataset iris.
#*# The data mining process in tdmClassify calls randomForest as the prediction model.
#*# It is called opts$NRUN=1 time with one random train-validation set splits.
#*# Therefore data frame res$allEval has one row
#*#
opts=tdmOptsDefaultsSet()                       # set all defaults for data mining process
gdObj <- tdmGraAndLogInitialize(opts);          # init graphics and log file

data(iris)
response.variables="Species"                    # names, not data (!)
input.variables=setdiff(names(iris),"Species")
opts$NRUN=1

idx_train = sample(nrow(iris))[1:110]
d_train=iris[idx_train,]
d_vali=iris[-idx_train,]
d_dis=iris[numeric(0),]
res <- tdmClassify(d_train,d_vali,d_dis,NULL,response.variables,input.variables,opts)

cat("\n")
print(res$allEVAL)
