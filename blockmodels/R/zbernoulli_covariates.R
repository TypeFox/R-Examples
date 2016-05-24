sigmo <- function(x) {1/(1+exp(-x))}

BM_bernoulli_covariates <- setRefClass("BM_bernoulli_covariates",
    contains = "scalar_model_with_covariates",
    methods = list(
        initialize = function(membership_type,adj,covariates,...)
        {
            if(class(covariates)=="matrix")
            {
                covariates_list <- list(covariates)
            }
            else
            {
                covariates_list <- covariates
            }
            .self$initFields(membership_name = membership_type,
                            model_name = "bernoulli_covariates",
                            adj = adj,
                            covariates = covariates_list,
                            ...)
            .self$postinit()
        },
        postinit = function()
        {
            callSuper()
        },
        plot_parameters = function(Q)
        {
            matrixplot(.self$plot_transform(model_parameters[[Q]]$m))
        },
        prediction = function(Q)
        {
            B <- matrix(0,nrow(adj),ncol(adj))
            for(k in 1:length(covariates))
            {
                B <- B + model_parameters[[Q]]$beta[k] * covariates[[k]]
            }

            if(membership_name=='LBM')
            {
                return(
                    sigmo(
                        memberships[[Q]]$Z1
                        %*%
                        model_parameters[[Q]]$m
                        %*%
                        t(memberships[[Q]]$Z2)
                        + B
                    )
                )
            }
            else
            {
                return(
                    sigmo(
                        memberships[[Q]]$Z
                        %*%
                        model_parameters[[Q]]$m
                        %*%
                        t(memberships[[Q]]$Z)
                        + B
                    )
                )
            }
        }
    )
)

