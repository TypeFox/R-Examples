Viterbi <- function (object, ...){
    if (class(object)=="numeric")
        do.call(Viterbihmm, list(object, ...))
    else if (class(object)=="mmpp") stop("Viterbi does not yet have a method for objects of class 'mmpp'.")
    else UseMethod("Viterbi")
}

