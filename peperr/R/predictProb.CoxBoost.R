predictProb.CoxBoost <- function (object, response, x, times, complexity, ...) 
{
    if (is.list(complexity)) {
        predict(object, type = "risk", newdata = x, times = times, 
            at.step = complexity$stepno)
    }
    else {
        predict(object, type = "risk", newdata = x, times = times, 
            at.step = complexity)
    }
}




