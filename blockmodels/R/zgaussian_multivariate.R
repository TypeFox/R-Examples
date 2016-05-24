
BM_gaussian_multivariate <- setRefClass("BM_gaussian_multivariate",
    contains = "multivariate_model",
    methods = list(
        initialize = function(membership_type,adj,...)
        {
            .self$initFields(membership_name = membership_type,
                            model_name = "gaussian_multivariate",
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
            iL <- chol(solve(model_parameters[[Q]]$Sigma))
            pred <- .self$prediction(Q)
            res <- list()
            for(k in 1:length(adj))
            {
                res[[k]] <- (adj[[k]]-pred[[k]])
            }

            normalized_residual <- list()
            for(k in 1:length(adj))
            {
                normalized_residual[[k]]<-matrix(0,nrow(adj[[1]]),ncol(adj[[1]]))
            }

            for(k1 in 1:length(adj))
            {
                for(k2 in 1:length(adj))
                {
                    normalized_residual[[k1]] <- normalized_residual[[k1]] + iL[k1,k2] * res[[k2]]
                }
            }
            return(normalized_residual)
        } 
    )
)

