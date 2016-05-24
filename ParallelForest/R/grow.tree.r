#------------------------------------------------------------------------------
#   Defines an R function to grow a tree with the compiled underlying
#       Fortran base. Returns an object of class tree.
#   Copyright (C) 2014  Bertram Ieong
#   No warranty provided.
#------------------------------------------------------------------------------


grow.tree = function(formula, data, subset, weights, na.action,
    impurity.function = "gini", model = FALSE, x = FALSE, y = FALSE,
    min_node_obs, max_depth){


    ### Input Assertions ###
    if(length(min_node_obs)!=1) stop ("min_node_obs must be a scalar.")
    if(length(max_depth)!=1)    stop ("max_depth must be a scalar.")




    if(min_node_obs<1) stop ("min_node_obs must be at least 1.")
    if(max_depth<0)    stop ("max_depth must be at least 0.")




    if(impurity.function!="gini"){
        stop("Only the Gini impurity function is currently supported.")
    }


    ### Create design matrix and dependent variable vector ###
    # create model frame #
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

    ytrain.tof = as.integer(ytrain)
    xtrain.tof = as.matrix(xtrain)
    storage.mode(xtrain.tof) = "double"

    # assert that Y must be 0 or 1
    y.unique.sorted = sort(unique(ytrain.tof))
    if(length(y.unique.sorted)<2) stop("Dependent variable must have two classes.")
    if(length(y.unique.sorted)>2) stop(paste("Dependent variable can only have two classes.",
        "Support for more classes may be implemented in a future version of this package"))

    if(sum(y.unique.sorted==c(0,1))!=2) stop(paste("Dependent variable must be automatically", 
        "coercible to classes 0 and 1. Please refactor the dependent variable to 0 and 1.",
        "More flexible automatic coercion may be implemented in a future version of this package."))

    # get data size
    n = nrow(xtrain)
    p = ncol(xtrain)


    ### Fit tree with Fortran compiled program ###

    # determine the maximum possible number of nodes with the given max depth for the 
    # fitted tree, which determines the length of the padded array that the Fortran
    # subroutine should return
    TOP_NODE_NUM = 0
    retlen = 2^(max_depth + 1 - TOP_NODE_NUM) - 1

    # check feasibility of passing Fortran to R results through memroy
    if((as.integer(retlen) > .Machine$integer.max) | (is.na(as.integer(retlen)))){
        stop(paste("grow.tree currently does not support",
            "inputs where (2^(max_depth + 1) - 1) exceeds",
            .Machine$integer.max,
            ". Support for larger values will be added in the next version of this package."))
    }


    # send to Fortran wrapper to grow forest
    ret = .Fortran("grow_tree_wrapper",
        n=as.integer(n), p=as.integer(p),
        xtrain=xtrain.tof, ytrain=ytrain.tof,
        min_node_obs=as.integer(min_node_obs), max_depth=as.integer(max_depth), 
        retlen=as.integer(retlen),
        tag_padded=integer(retlen),
        tagparent_padded=integer(retlen),
        tagleft_padded=integer(retlen),
        tagright_padded=integer(retlen),
        is_topnode_padded=integer(retlen),
        depth_padded=integer(retlen),
        majority_padded=integer(retlen),
        has_subnodes_padded=integer(retlen),
        splitvarnum_padded=integer(retlen),
        splitvalue_padded=double(retlen),
        numnodes=integer(1)
        )

    # unpad returned arrays and put everything into a tree object
    flattened.nodes = data.frame(
        tag=ret$tag_padded[1:ret$numnodes],
        tagparent=ret$tagparent_padded[1:ret$numnodes],
        tagleft=ret$tagleft_padded[1:ret$numnodes],
        tagright=ret$tagright_padded[1:ret$numnodes],
        is_topnode=ret$is_topnode_padded[1:ret$numnodes],
        depth=ret$depth_padded[1:ret$numnodes],
        majority=ret$majority_padded[1:ret$numnodes],
        has_subnodes=ret$has_subnodes_padded[1:ret$numnodes],
        splitvarnum=ret$splitvarnum_padded[1:ret$numnodes],
        splitvalue=ret$splitvalue_padded[1:ret$numnodes]
        )

    fitted.tree = new("tree",
        n=ret$n, p=ret$p,
        min_node_obs=ret$min_node_obs, max_depth=ret$max_depth,
        numnodes=ret$numnodes,
        flattened.nodes=flattened.nodes,
        fmla=formula
        )

    return(fitted.tree)
}

