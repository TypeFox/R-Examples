findUL <- function (step = 0.01, initStep = 0, fun, MLE, level = 3.84, 
    ...) 
{
    value <- 0
    step1 <- step
    Lbeta <- MLE - initStep
    for (i in 1:8) {
        while (value < level) {
            Lbeta <- Lbeta - step1
            value <- fun(Lbeta, ...)$"-2LLR"
        }
        Lbeta <- Lbeta + step1
        step1 <- step1/10
        value <- fun(Lbeta, ...)$"-2LLR"
    }
    value1 <- value
    value <- 0
    Ubeta <- MLE + initStep
    for (i in 1:8) {
        while (value < level) {
            Ubeta <- Ubeta + step
            value <- fun(Ubeta, ...)$"-2LLR"
        }
        Ubeta <- Ubeta - step
        step <- step/10
        value <- fun(Ubeta, ...)$"-2LLR"
    }
	if ( (value1>level)|(value>level) ) warning("Something wrong. Check the MLE and step inputs.")
    return(list(Low = Lbeta, Up = Ubeta, FstepL = step1, FstepU = step, 
        Lvalue = value1, Uvalue = value))
}