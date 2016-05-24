## ----eval=FALSE----------------------------------------------------------
#  install.packages("wsrf")

## ----eval=FALSE----------------------------------------------------------
#  install.packages("wsrf",
#                   type="source",
#  				 configure.args="--enable-c11=yes")

## ----eval=FALSE----------------------------------------------------------
#  install.packages("wsrf",
#                   configure.args="--enable-c11=no")

## ----eval=FALSE----------------------------------------------------------
#  install.packages("wsrf",
#                   type="source",
#                   configure.args="--with-boost-include=<Boost include path>
#                                   --with-boost-lib=<Boost lib path>")

## ----usage_load, message=FALSE-------------------------------------------
library("rattle")
ds <- weather
dim(ds)
names(ds)

## ----usage_prepare-------------------------------------------------------
target <- "RainTomorrow"
ignore <- c("Date", "Location", "RISK_MM")
(vars <- setdiff(names(ds), ignore))
dim(ds[vars])

## ----message=FALSE-------------------------------------------------------
library("randomForest")
if (sum(is.na(ds[vars]))) ds[vars] <- na.roughfix(ds[vars])
ds[target] <- as.factor(ds[[target]])
(tt <- table(ds[target]))

## ------------------------------------------------------------------------
(form <- as.formula(paste(target, "~ .")))

## ------------------------------------------------------------------------
seed <- 42
set.seed(seed)
length(train <- sample(nrow(ds), 0.7*nrow(ds)))
length(test <- setdiff(seq_len(nrow(ds)), train))

## ----eval=FALSE----------------------------------------------------------
#  wsrf(formula,
#       data,
#       ntrees=500,
#       nvars=NULL,
#       weights=TRUE,
#       parallel=TRUE)

## ----usage_build_by_default, message=FALSE-------------------------------
library("wsrf")
model.wsrf.1 <- wsrf(form, data=ds[train, vars])
print(model.wsrf.1)
print(model.wsrf.1, 1)  # Print tree 1.

## ----usage_evaluate------------------------------------------------------
cl <- predict(model.wsrf.1, newdata=ds[test, vars], type="class")
actual <- ds[test, target]
(accuracy.wsrf <- sum(cl == actual, na.rm=TRUE)/length(actual))

## ----usage_build_another, message=FALSE----------------------------------
# Here we build another model without weighting.
model.wsrf.2 <- wsrf(form, data=ds[train, vars], weights=FALSE)
print(model.wsrf.2)

## ----usage_subset_combine------------------------------------------------
submodel.wsrf <- subset.wsrf(model.wsrf.1, 1:150)
print(submodel.wsrf)
bigmodel.wsrf <- combine.wsrf(model.wsrf.1, model.wsrf.2)
print(bigmodel.wsrf)

## ----usage_build_on_cluster, eval=FALSE----------------------------------
#  servers <- paste0("node", 31:40)
#  model.wsrf.3 <- wsrf(form, data=ds[train, vars], parallel=servers)

