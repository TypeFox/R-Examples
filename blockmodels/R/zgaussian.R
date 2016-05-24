
BM_gaussian <- setRefClass("BM_gaussian",
    contains = "scalar_model",
    methods = list(
        initialize = function(membership_type,adj,...)
        {
            .self$initFields(membership_name = membership_type,
                            model_name = "gaussian",
                            adj = adj,
                            ...)
            .self$postinit()
        },
        plot_parameters = function(Q)
        {
            matrixplot(model_parameters[[Q]]$mu)
        },
        prediction = function(Q)
        {
            if(membership_name=='LBM')
            {
                return(
                    memberships[[Q]]$Z1
                     %*%
                    model_parameters[[Q]]$mu
                     %*%
                    t(memberships[[Q]]$Z2)
                )
            }
            else
            {
                return(
                    memberships[[Q]]$Z
                     %*%
                    model_parameters[[Q]]$mu
                     %*%
                    t(memberships[[Q]]$Z)
                )
            }
        }
    )
)

