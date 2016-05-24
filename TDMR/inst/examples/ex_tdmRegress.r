#*# This example shows a simple data mining process (phase 1 of TDMR) for regression on
#*# dataset iris.
#*# The data mining process in tdmRegress calls randomForest as the prediction model.
#*# It is called  for 2 response variables. Therefore, the data frames allRMAE and allRMSE 
#*# have 2 rows.
#*#
opts=tdmOptsDefaultsSet()                       # set all defaults for data mining process
gdObj <- tdmGraAndLogInitialize(opts);          # init graphics and log file

data(iris)
response.variables=c("Petal.Length","Petal.Width")                # names, not data (!)
input.variables=setdiff(names(iris),response.variables)
opts$rgain.type="rmae"
opts$NRUN=1

idx_train = sample(nrow(iris))[1:110]
d_train=iris[idx_train,]
d_vali=iris[-idx_train,]
res <- tdmRegress(d_train,d_vali,NULL,response.variables,input.variables,opts)

print(res$allRMAE)
print(res$allRMSE)
