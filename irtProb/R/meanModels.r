`meanModels` <-
function(modelShow, statistics=c("T","S")) {
 res           <- data.frame(MODEL = modelShow$MODEL, modelShow[,statistics])
 res           <- aggregate(res[,-1], list(MODEL=res$MODEL), mean, na.rm=TRUE)
 rownames(res) <- res$MODEL
 return(res[,-1])
 }


