GLMBoost <- function(x,y,penalty=length(y),standardize=TRUE,...) {
    res <- GAMBoost(x=NULL,x.linear=x,y=y,penalty.linear=penalty,standardize.linear=standardize,...)
    class(res) <- c("GLMBoost",class(res))
    res
}

predict.GLMBoost <- function(object,newdata=NULL,...) {
    predict.GAMBoost(object,newdata=NULL,newdata.linear=newdata,...)
}

cv.GLMBoost <- function(x,y,penalty=length(y),just.criterion=TRUE,...) {
    res <- cv.GAMBoost(x=NULL,x.linear=x,y=y,penalty.linear=penalty,just.criterion=just.criterion,...)
    if (!just.criterion) class(res) <- c("GLMBoost",class(res))
    res
}

optimGLMBoostPenalty <- function(x,y,start.penalty=length(y),just.penalty=FALSE,...) {
    res <- optimGAMBoostPenalty(x=NULL,x.linear=x,y=y,start.penalty=start.penalty,
                                just.penalty=just.penalty,which.penalty="linear",...)
    if (!just.penalty) class(res) <- c("GLMBoost",class(res))
    res
}