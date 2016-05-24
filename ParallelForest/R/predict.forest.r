#------------------------------------------------------------------------------
#   Defines a predict method for the class forest to make predictions on
#       fitted forests.
#   Copyright (C) 2014  Bertram Ieong
#   No warranty provided.
#------------------------------------------------------------------------------

# method to be used for predict function in forest S4 class
predict.forest = function(object, newdata, ...){

    ### Create design matrix ###
    fmla.str = deparse(object@fmla)
    Yvar.str = strsplit(fmla.str," ~ ")[[1]][1]

    newdata.tmp = newdata
    newdata.tmp[,Yvar.str] = rep(-1,nrow(newdata.tmp))

    xtest = model.frame(object@fmla, data=newdata.tmp)[,-1]
    xtest.tof = as.matrix(xtest)
    storage.mode(xtest.tof) = "double"

    # get new data size #
    n.new = as.integer(nrow(xtest))
    p = as.integer(ncol(xtest))

    if(p != object@p){
        stop("New data has different number of variables than training data.")
    }

    ### Send to compiled Fortran wrapper to get forest prediction ###
    retpred = .Fortran("predict_forest_wrapper",
        (object@flattened.nodes)$tag,
        (object@flattened.nodes)$tagparent,
        (object@flattened.nodes)$tagleft,
        (object@flattened.nodes)$tagright,
        (object@flattened.nodes)$is_topnode,
        (object@flattened.nodes)$depth,
        (object@flattened.nodes)$majority,
        (object@flattened.nodes)$has_subnodes,
        (object@flattened.nodes)$splitvarnum,
        (object@flattened.nodes)$splitvalue,
        object@numnodes,
        object@numboots,
        n.new,
        p,
        xtest.tof,
        ynew_pred_01=integer(n.new)
        )

    ### Return vector with forest prediction ###

    ynew_pred = do.call(prep.depvar.out, append(list(retpred$ynew_pred_01), object@depvar.restore.info))
    return(ynew_pred)
}

# set the predict function with the method defined above for the forest S4 class
setMethod(f = "predict", signature = "forest", definition = predict.forest)
