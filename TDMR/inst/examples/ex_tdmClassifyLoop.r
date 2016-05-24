#*# --------- demo/demo00-0classif.r ---------
#*# This demo shows a simple data mining process (phase 1 of TDMR) for classification on
#*# dataset iris.
#*# The data mining process in tdmClassifyLoop calls randomForest as the prediction model.
#*# It is called opts$NRUN=2 times with different random train-validation set splits.
#*# Therefore data frame result$Err has two rows
#*#
opts=tdmOptsDefaultsSet()                       # set all defaults for data mining process
gdObj <- tdmGraAndLogInitialize(opts);          # init graphics and log file

data(iris)
response.variables="Species"                    # names, not data (!)
input.variables=setdiff(names(iris),"Species")

result = tdmClassifyLoop(iris,response.variables,input.variables,opts)

print(result$Err)
