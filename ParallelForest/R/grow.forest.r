#------------------------------------------------------------------------------
#   Defines an R function to grow a forest with the compiled underlying
#       Fortran base. Returns an object of class forest.
#   Copyright (C) 2014  Bertram Ieong
#   No warranty provided.
#------------------------------------------------------------------------------


grow.forest = function(formula, data, subset, na.action,
    impurity.function = "gini", model = FALSE, x = FALSE, y = FALSE,
    min_node_obs, max_depth, 
    numsamps, numvars, numboots=5){


    ### Create design matrix and dependent variable vector ###
    # create model frame #
    if(missing(data)) data=environment(formula)

    if(missing(subset) & missing(na.action)){
        m = model.frame(formula, data=data)
    } else if(missing(subset) & !missing(na.action)){
        m = model.frame(formula, data=data, na.action=na.action)
    } else if(!missing(subset) & missing(na.action)){
        m = model.frame(formula, data=data, subset=subset)
    } else if(!missing(subset) & !missing(na.action)){
        m = model.frame(formula, data=data, subset=subset, na.action=na.action)
    } else {
        stop("Error.")
    }

    # create matrices to be fed to Fortran compiled program
    ytrain = m[,1]
    xtrain = m[,-1]

    retlist = prep.depvar.in(ytrain)

    myvars = names(retlist) %in% "y.prepped"
    depvar.restore.info = retlist[!myvars]

    ytrain.tof = retlist$y.prepped
    xtrain.tof = as.matrix(xtrain)
    storage.mode(xtrain.tof) = "double"


    # get data size
    n = nrow(xtrain)
    p = ncol(xtrain)


    ### Put in default values if missing ###
    if(missing(min_node_obs)){
        min_node_obs = ceiling(n/200)
    }
    if(missing(max_depth)){
        max_depth = ceiling(n/20000)
    }
    if(missing(numsamps)){
        numsamps = ceiling(n/2)
    }
    if(missing(numvars)){
        numvars = ceiling(0.70 * p)
    }


    ### Input Assertions ###
    if(length(min_node_obs)!=1) stop ("min_node_obs must be a scalar.")
    if(length(max_depth)!=1)    stop ("max_depth must be a scalar.")
    if(length(numsamps)!=1)     stop ("numsamps must be a scalar.")
    if(length(numvars)!=1)      stop ("numvars must be a scalar.")
    if(length(numboots)!=1)     stop ("numboots must be a scalar.")

    if(min_node_obs<1) stop ("min_node_obs must be at least 1.")
    if(max_depth<0)    stop ("max_depth must be at least 0.")
    if(numsamps<1)     stop ("numsamps must be at least 1.")
    if(numvars<1)      stop ("numvars must be at least 1.")
    if(numboots<1)     stop ("numboots must be at least 1.")

    if(impurity.function!="gini"){
        stop("Only the Gini impurity function is currently supported.")
    }


    ### Fit forest with Fortran compiled program ###

    # determine the maximum possible number of nodes with the given max depth for the 
    # fitted tree, which determines the length of the padded array that the Fortran
    # subroutine should return
    TOP_NODE_NUM = 0
    retlen = 2^(max_depth + 1 - TOP_NODE_NUM) - 1

    # check feasibility of passing Fortran to R results through memory
    if((as.integer(retlen*numboots) > .Machine$integer.max) | (is.na(as.integer(retlen*numboots)))){
        stop(paste("grow.forest currently does not support",
            "inputs where (2^(max_depth + 1) - 1) * numboots exceeds",
            .Machine$integer.max,
            ". Support for larger values will be added in the next version of this package."))
    }


    # send to Fortran wrapper to grow forest
    ret = .Fortran("grow_forest_wrapper",
        n=as.integer(n), p=as.integer(p),
        xtrain=xtrain.tof, ytrain=ytrain.tof,
        min_node_obs=as.integer(min_node_obs), max_depth=as.integer(max_depth), 
        retlen=as.integer(retlen),
        numsamps=as.integer(numsamps),
        numvars=as.integer(numvars),
        numboots=as.integer(numboots),
        tag_padded=integer(retlen*numboots),
        tagparent_padded=integer(retlen*numboots),
        tagleft_padded=integer(retlen*numboots),
        tagright_padded=integer(retlen*numboots),
        is_topnode_padded=integer(retlen*numboots),
        depth_padded=integer(retlen*numboots),
        majority_padded=integer(retlen*numboots),
        has_subnodes_padded=integer(retlen*numboots),
        splitvarnum_padded=integer(retlen*numboots),
        splitvalue_padded=double(retlen*numboots),
        numnodes=integer(numboots)
        )

    # unpad returned arrays and put everything into a forest object
    flattened.nodes = data.frame(
        tag=ret$tag_padded[1:sum(ret$numnodes)],
        tagparent=ret$tagparent_padded[1:sum(ret$numnodes)],
        tagleft=ret$tagleft_padded[1:sum(ret$numnodes)],
        tagright=ret$tagright_padded[1:sum(ret$numnodes)],
        is_topnode=ret$is_topnode_padded[1:sum(ret$numnodes)],
        depth=ret$depth_padded[1:sum(ret$numnodes)],
        majority=ret$majority_padded[1:sum(ret$numnodes)],
        has_subnodes=ret$has_subnodes_padded[1:sum(ret$numnodes)],
        splitvarnum=ret$splitvarnum_padded[1:sum(ret$numnodes)],
        splitvalue=ret$splitvalue_padded[1:sum(ret$numnodes)]
        )

    fitted.forest = new("forest",
        n=ret$n, p=ret$p,
        min_node_obs=ret$min_node_obs, max_depth=ret$max_depth,
        numsamps=ret$numsamps,
        numvars=ret$numvars,
        numboots=ret$numboots,
        numnodes=ret$numnodes,
        flattened.nodes=flattened.nodes,
        fmla=formula,
        depvar.restore.info=depvar.restore.info
        )

    ### Store model frame, x, and y if requested ###
    if(model) fitted.forest@model = m
    if(x) fitted.forest@x     = xtrain
    if(y) fitted.forest@y     = ytrain

    ### Return fitted forest object ###
    return(fitted.forest)
}

