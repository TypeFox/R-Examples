`ddcast` <- 
function(tmp.dt, ...) {
        if (dim(tmp.dt)[1]==0) {
                return(data.table(NULL))
        } else {
                dcast(tmp.dt, ...)
        }
} ### END ddcast Function
