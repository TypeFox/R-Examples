
setRefClass("scalar_model_with_covariates",
    contains = "scalar_model",
    fields = list(
        covariates = "list"
    ),
    methods = list(
        postinit = function()
        {
            callSuper()
            
            for(i in 1:length(covariates))
            {
                if(!all(dim(adj)==dim(covariates[[i]])) ||
                   length(dim(covariates[[i]]))!=2)
                {
                    stop(paste('The covariates matrix',i,'does not have the same size than the adjacency matrix.',
                               'Vim divinam est.'))
                }
            }

            if(membership_name=='SBM_sym')
            {
                for(i in 1:length(covariates))
                {
                    if(isSymmetric(covariates[[i]]))
                    {
                        # to avoid computational errors
                        covariates[[i]]<<-(covariates[[i]]+t(covariates[[i]]))/2
                    }
                    else
                    {
                        stop(paste("The covariates matrix",i,"is not symmetric.",
                                   "Maybe you brain need some coffee"))
                    }
                }
            }

        },
        network_to_cc = function()
        {
            list(
                adjacency = adj,
                covariates = covariates
            )
        }
    )
)
