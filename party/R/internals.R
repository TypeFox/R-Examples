
party_intern <- function(..., fun = c("R_TreeGrow", "R_get_nodeID",
    "R_getpredictions", "R_modify_response", "R_remove_weights", 
    "ctreedpp", "newinputs")) {

    fun <- match.arg(fun)
    do.call(fun, list(...))
}
