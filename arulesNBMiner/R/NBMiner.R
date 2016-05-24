NBMiner <- function(data, parameter, control = NULL) {

    data <- as(data, "transactions")
    #parameter <- as(parameter, "NBMinerParameter")
    control <- as(control, "NBMinerControl")

    if(control@verbose) {
        ## print parameter
        cat("\nparameter specification:\n")
        print(parameter)
        cat("\nalgorithmic control:\n")
        print(control)
        cat("\n")
    }

    ## create DB
    db <- .jnew("SparseSetOfItemsets", data@data@i, data@data@p, 
        dim(data)[2])
    
    ## call NBMiner
    miner <- .jnew("NBMiner")
    result <- .jcall(miner, "LR_result;", "R_mine", db,
            parameter@pi, parameter@theta, parameter@a, parameter@k, 
            parameter@n, parameter@maxlen, 
            parameter@rules,
            control@verbose, control@debug)

    ## get result
    .as_itemMatrix <- function(x) {
        m <- new("ngCMatrix", 
            i = .jcall(x, "[I", "getI"),
            p = .jcall(x, "[I", "getP"),
            Dim = c(.jcall(x, "I", "getItems"), 
                    .jcall(x, "I", "size"))
        )

        new("itemMatrix", data = m, itemInfo = itemInfo(data))
    }

    ## get precision
    precision <- .jcall(result, "[D", "getPrecision")
    
    ## encode as rules/itemsets
    res <- if(parameter@rules) new("rules", 
        lhs = .as_itemMatrix(.jcall(result, "LSparseSetOfItemsets;", "getLhs")),
        rhs = .as_itemMatrix(.jcall(result, "LSparseSetOfItemsets;", "getRhs")),
        quality = data.frame(precision=precision))
    else new("itemsets", 
        items = .as_itemMatrix(.jcall(result, 
                "LSparseSetOfItemsets;", "getItems")),
        quality = data.frame(precision=precision))

    ## remove itemsets/rules that are too short
    if(parameter@minlen>1) res <- res[size(res) >= parameter@minlen]
    res
}
