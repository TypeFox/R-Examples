print.extraTrees <- function(x, ...) {
    cat( "ExtraTrees:\n" )
    cat( sprintf(" - # of trees: %d\n", x$ntree) )
    cat( sprintf(" - node size:  %d\n", x$nodesize) )
    cat( sprintf(" - # of dim:   %d\n", x$ndim) )
    cat( sprintf(" - # of tries: %d\n", x$mtry) )
    if (x$factor) {
    	type="factor (classification)"
    } else {
    	type="numeric (regression)"
    }
    cat( sprintf(" - type:       %s\n", type) )
    cat( sprintf(" - multi-task: %s\n", ifelse(x$multitask, "yes", "no")) )
}

