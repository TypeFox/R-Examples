mtest.gv <-
function(gv) {
    # test if 
    if (class(gv) != 'gv') stop("Object must be of class 'gv'")

    if (length(gv$model) == 1) {
        if (is.na(gv$model)) {
            mtest <- FALSE
        } else {
            stop("Model in 'gv' has no parameters")
        }
    } else {
        if (gv$model$type %in% 1:3 & length(gv$model) == 4) {
            mtest <- TRUE
        } else if (gv$model$type == 4 & length(gv$model) == 5) {
            mtest <- TRUE
        }
    }
    mtest
}
