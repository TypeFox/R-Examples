`export` <-
function(x, file = "", ...) {
    
    export.args <- list(...)
    Call <- match.call(expand.dots = TRUE)
    # Callist <- as.list(Call)
    # as.list(Call) e acelasi lucru cu list(...)
    
    caseid <- "cases"
    if (any(names(export.args) == "caseid")) {
        caseid <- export.args[["caseid"]]
        Call[["caseid"]] <- NULL
    }
    
    
    if (!missing(x)) {
        if (is.data.frame(x)) {
            if (any(rownames(x) != seq(nrow(x)))) {
                if (all(colnames(x) != caseid)) {
                    x <- cbind("cases" = rownames(x), x)
                    names(x)[1] <- caseid
                }
            }
        }
    }
    
    if (any(names(export.args) == "row.names")) {
        warning("The argument \"row.names\" is set to FALSE by default.", domain = NA)
    }
    
    if (any(names(export.args) == "sep")) {
        if (export.args[["sep"]] == "tab") {
            export.args[["sep"]] <- "\t"
        }
        Call[["sep"]] <- export.args[["sep"]]
    }
    else {
        Call[["sep"]] <- ","
    }
    
    if (any(names(export.args) == "col.names")) {
        Call[["col.names"]] <- export.args[["col.names"]]
    }
    
                                                       
    Call[[1L]] <- as.name("write.table")
    Call[[2L]] <- as.name("x")
    Call[[3L]] <- file
    Call[["row.names"]] <- FALSE
    
    eval(Call)
}
