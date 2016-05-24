
BM_gaussian_multivariate_independent_homoscedastic <- setRefClass("BM_gaussian_multivariate_independent_homoscedastic",
    contains = "multivariate_model",
    methods = list(
        initialize = function(membership_type,adj,...)
        {
            .self$initFields(membership_name = membership_type,
                            model_name = "gaussian_multivariate_independent_homoscedastic",
                            adj = adj,
                            ...)
            .self$postinit()
        },
        plot_parameters = function(Q)
        {
        },
        prediction = function(Q)
        {
            if(membership_name=='LBM')
            {
                pred <- list()
                for(k in 1:length(adj))
                {
                    pred[[k]]<- memberships[[Q]]$Z1 %*% model_parameters[[Q]]$mu[[k]] %*% t(memberships[[Q]]$Z2)
                }
                return(pred)
            }
            else
            {
                pred <- list()
                for(k in 1:length(adj))
                {
                    pred[[k]] <- memberships[[Q]]$Z %*% model_parameters[[Q]]$mu[[k]] %*% t(memberships[[Q]]$Z)
                }
                return(pred)
            }
        },
        residual = function(Q)
        {
            pred <- .self$prediction(Q)
            ret <- list()
            for(k in 1:length(adj))
            {
                ret[[k]] <- adj[[k]]-pred[[k]]
            }
            return(ret)
        } 
    )
)

