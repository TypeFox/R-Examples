suppressMessages(library(wsrf))
suppressMessages(library(rattle))
suppressMessages(library(randomForest))

# prepare parameters
library("wsrf")
library("rattle")
library("randomForest")
ds <- weather
target <- "RainTomorrow"
ignore <- c("Date", "Location", "RISK_MM")
vars <- setdiff(names(ds), ignore)
if (sum(is.na(ds[vars]))) ds[vars] <- na.roughfix(ds[vars])
ds[target] <- as.factor(ds[[target]])
(form <- as.formula(paste(target, "~ .")))
set.seed(500)
length(train <- sample(nrow(ds), 0.7*nrow(ds)))
length(test  <- setdiff(seq_len(nrow(ds)), train))


# build model
model.wsrf       <- wsrf(form, data=ds[train, vars])
model.wsrf.nw    <- wsrf(form, data=ds[train, vars], weights=FALSE)
model.wsrf.nw.vi <- wsrf(form, data=ds[train, vars], weights=FALSE, importance=TRUE)
model.subset     <- subset.wsrf(model.wsrf, 1:200)
model.combine    <- combine.wsrf(model.wsrf, model.wsrf.nw)


# evaluate

# Note:
# 32bit system and 64bit system will have different results, however,
# if random seed is fixed, the same results will be presented in the
# same system.

cl         <- predict(model.wsrf,     newdata=ds[test, vars], type="class")
cl.nw      <- predict(model.wsrf.nw,  newdata=ds[test, vars], type="class")
cl.subset  <- predict(model.subset,   newdata=ds[test, vars], type="class")
cl.combine <- predict(model.combine,  newdata=ds[test, vars], type="class")
