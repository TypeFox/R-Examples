
BM_poisson_covariates <- setRefClass("BM_poisson_covariates",
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
                            model_name = "poisson_covariates",
                            adj = adj,
                            covariates = covariates_list,
                            ...)
            .self$postinit()
        },
        plot_transform = function(x) {log(1+x)},
        plot_parameters = function(Q)
        {
            matrixplot(.self$plot_transform(model_parameters[[Q]]$lambda))
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
                    (
                        memberships[[Q]]$Z1
                        %*%
                        model_parameters[[Q]]$lambda
                        %*%
                        t(memberships[[Q]]$Z2)
                    ) * exp(B)
                )
            }
            else
            {
                return(
                    (
                        memberships[[Q]]$Z
                        %*%
                        model_parameters[[Q]]$lambda
                        %*%
                        t(memberships[[Q]]$Z)
                    ) * exp(B)
                )
            }
        }
    )
)

