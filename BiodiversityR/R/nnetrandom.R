`nnetrandom` <-
function(formula,data,tries=10,leave.one.out=F,...){
    nnet.1 <- function(formula,data,tries,...) {
        optimal <- Inf
        values <- numeric(length=tries)
        for (i in 1:tries) {
            temp.result <- nnet::nnet(formula=formula,...)
            values[i] <- temp.result$value
            if (temp.result$value < optimal) {
                final.result <- temp.result
                optimal <- temp.result$value
            }
        }
        final.result$range <- summary(values)
        return(final.result)
    }
    result <- nnet.1(formula=formula,data=data,tries=tries,...)
    result$tries <- tries
    if (leave.one.out==T) {
        predictresult <- character(length=nrow(data))
        respvar <- all.vars(formula)[1]
        for (i in 1:nrow(data)) {
            data1 <- data[-i,]
            result1 <- nnet.1(formula=formula,data=data1,tries=tries,...)
            predictresult[i] <- predict(result1,data=data,type="class")[i]
        }
        result$CV <- predictresult
        result$successful <- predictresult==data[,respvar]
    }
    return(result)
}
