# Display clustering results
setMethod("show", signature(object="APResult"),
    function(object)
    {
        cat("\nAPResult object\n")

        if (!is.finite(object@l) || !is.finite(object@it))
            stop("object is not result of an affinity propagation run; ",
                 "it is pointless to create 'APResult' objects yourself.")

        cat("\nNumber of samples     = ", object@l, "\n")
        if (length(object@sel) > 0)
        {
            cat("Number of sel samples = ", length(object@sel),
                paste("   (", round(100*length(object@sel)/object@l,1),
                      "%)\n", sep=""))
            cat("Number of sweeps      = ", object@sweeps, "\n")
        }
        cat("Number of iterations  = ", object@it, "\n")
        cat("Input preference      = ", object@p, "\n")
        cat("Sum of similarities   = ", object@dpsim, "\n")
        cat("Sum of preferences    = ", object@expref, "\n")
        cat("Net similarity        = ", object@netsim, "\n")
        cat("Number of clusters    = ", length(object@exemplars), "\n\n")

        if (length(object@exemplars) > 0)
        {
            if (length(names(object@exemplars)) == 0)
            {
                cat("Exemplars:\n")
                cat(object@exemplars, fill=TRUE, labels="  ")
                cat("Clusters:\n")

                for (i in 1:length(object@exemplars))
                {
                    cat("   Cluster ", i, ", exemplar ",
                        object@exemplars[i], ":\n", sep="")
                    cat(object@clusters[[i]], fill=TRUE, labels="     ")
                }
            }
            else
            {
                cat("Exemplars:\n")
                cat(names(object@exemplars), fill=TRUE, labels="  ")
                cat("Clusters:\n")

                for (i in 1:length(object@exemplars))
                {
                    cat("   Cluster ", i, ", exemplar ",
                        names(object@exemplars[i]), ":\n", sep="")
                    cat(names(object@clusters[[i]]), fill=TRUE, labels="     ")
                }
            }
        }
        else
        {
            cat("No clusters identified.\n")
        }
    }
)

setMethod("show", signature(object="ExClust"),
    function(object)
    {
        cat("\nExClust object\n")

        if (!is.finite(object@l))
            stop("object is not result of an exemplar-based clustering; ",
                 "it is pointless to create 'ExClust' objects yourself.")

        cat("\nNumber of samples   = ", object@l, "\n")
        cat("Number of clusters  = ", length(object@exemplars), "\n\n")

        if (length(object@exemplars) > 0)
        {
            if (length(names(object@exemplars)) == 0)
            {
                cat("Exemplars:\n")
                cat(object@exemplars, fill=TRUE, labels="  ")
                cat("Clusters:\n")

                for (i in 1:length(object@exemplars))
                {
                    cat("   Cluster ", i, ", exemplar ",
                        object@exemplars[i], ":\n", sep="")
                    cat(object@clusters[[i]], fill=TRUE, labels="     ")
                }
            }
            else
            {
                cat("Exemplars:\n")
                cat(names(object@exemplars), fill=TRUE, labels="  ")
                cat("Clusters:\n")

                for (i in 1:length(object@exemplars))
                {
                    cat("   Cluster ", i, ", exemplar ",
                        names(object@exemplars[i]), ":\n", sep="")
                    cat(names(object@clusters[[i]]), fill=TRUE, labels="     ")
                }
            }
        }
        else
        {
            cat("No clusters identified.\n")
        }
    }
)


setMethod("show", signature(object="AggExResult"),
    function(object)
    {
        cat("\nAggExResult object\n")

        if (!is.finite(object@l) || !is.finite(object@maxNoClusters))
            stop("object is not result of agglomerative clustering; ",
                 "it is pointless to create 'AggExResult' objects yourself.")

        cat("\nNumber of samples          = ", object@l, "\n")
        cat("Maximum number of clusters = ", object@maxNoClusters, "\n")
   }
)
