dist.multiPhylo <- function(x, method="geodesic", force.multi2di = FALSE,
        outgroup = NULL, convert.multifurcating = FALSE, use.random.resolution =
        FALSE, scale = NULL, verbose = FALSE)
{
    if(length(x) < 2)
        return(matrix())

    if(class(x) == "multiPhylo") # ideally, we will have this
    {
        # run checks if appropriate
        # separate it out into a vector of strings

        if(!is.null(outgroup))
        {
            x <- lapply(x, function(k)
                    {
                    if(class(k) == "phylo")
                    {
                        if(!is.rooted(k))
                            root(k, outgroup, resolve.root = TRUE)
                        else k
                    }
                    else NA
                    })
        }


        if(force.multi2di == TRUE)
        {
            x <- lapply(x, function(k) 
                  { 
                    if(class(k) == "phylo")
                    {
                    multi2di(k,random=use.random.resolution) 
                    }
                    else NA
                  }) 
        }
        else if(convert.multifurcating == TRUE) # won't resolve multifurcations at the root
        {
            x <- lapply(x, function(k) 
                  { 
                    if(class(k) == "phylo")
                    {
                        if(!is.binary.tree(k))
                            multi2di(k,random=use.random.resolution) 
                        else k 
                    }
                    else NA
                  }) 

        }

        if(!is.null(scale))
        {
            if(class(scale) == "phylo")
            {
                T <- sum(scale$edge.length)
                x <- lapply(x, function(k)
                        {
                            if(class(k) == "phylo")
                            {
                                k$edge.length = k$edge.length * (T /
                                    sum(k$edge.length))
                                k
                            }
                            else NA
                            })
            }
            else if(class(scale) == "numeric")
            {
                x <- lapply(x, function(k)
                        {
                            if(class(k) == "phylo")
                            {
                                k$edge.length = k$edge.length * (scale /
                                    sum(k$edge.length))
                                k
                            }
                            else NA
                            })
            }
            else
            {
                stop("Scale parameter not understood.\n")
            }
        }

        # do some sanity checks before we start out

        r <- lapply(x, function(k) 
                { 
                if(class(k) == "phylo")
                {
                    is.rooted(k) 
                }
                else NA
                })

        if(!all(as.logical(r), na.rm=TRUE))
        {
            stop("Some trees are not rooted. Specify an outgroup to fix this problem. All trees must be rooted.\n")
        }

        r <- lapply(x, function(k) { if(class(k) == "phylo") is.binary.tree(k) else NA })
        if(!all(as.logical(r), na.rm=TRUE))
        {
            stop("Some trees are not binary. All input trees must be strictly binary.\n")
        }

        # check to make sure all trees have the same tip labels
        tips = x[[1]]$tip.label
        r <- lapply(x, function(k) { if(class(k) == "phylo") setequal(k$tip.label,
                    tips) else NA})

        if(!all(as.logical(r), na.rm=TRUE))
        {
            stop("Not all trees have the same tips.\n")
        }

        # convert our list of class phylo to a list of strings
        treestrs <- lapply(x, function(k) { if(class(k) == "phylo")
                write.tree(k) else "" })
        
        method=tolower(method)

        method.id = pmatch(method, c("edgeset", "geodesic"))

        # call the C interface function and return the value automatically
        if(method.id == 1)
        {
            rv <- .Call("phycpp_bin_trees", treestrs, PACKAGE="distory")
        }
        else if(method.id == 2)
        {
            rv <- .Call("phycpp_compute_tree_distance_set", treestrs, as.logical(verbose), PACKAGE="distory")
        }
        else
        {
            stop("Specified method is not valid")
        }

        as.dist(rv)
    }

    else if(typeof(x) == "list")
    {
        if(class(x[[1]]) == "phylo") # a list of phylo's that for some reason is not classed as multiPhylo
        {
            class(x) <- "multiPhylo" # it already is basically a multiPhylo anyways - we'll mark it as such
            dist.multiPhylo(x, method=method, force.multi2di=force.multi2di, outgroup=outgroup,
                        convert.multifurcating=convert.multifurcating,
                        use.random.resolution=use.random.resolution,
                        scale=scale, verbose=verbose)
        }
        else if(class(x[[1]]) == "character") # a list of strings, presuming one tree each, properly terminated
        {
            # read with /ape/, run checks, dump back
            t <- paste(x, sep="", collapse="")
            k <- read.tree(text=t)
            dist.multiPhylo(x, method=method, force.multi2di=force.multi2di, outgroup=outgroup,
                        convert.multifurcating=convert.multifurcating,
                        use.random.resolution=use.random.resolution,
                        scale=scale, verbose=verbose)
        }
    }
    else if(class(x) == "character") # this is for one string containing multiple trees
    {
        # read with ape and dump back to a vector of strings
        k <- read.tree(text=x)

        # call this to process it properly
        dist.multiPhylo(x, method=method, force.multi2di=force.multi2di, outgroup=outgroup,
                        convert.multifurcating=convert.multifurcating,
                        use.random.resolution=use.random.resolution,
                        scale=scale, verbose=verbose)
    }
    else
    {
        stop("Cannot coerce the argument into a usable type.")
    }
}

