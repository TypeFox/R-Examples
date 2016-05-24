
BM_poisson <- setRefClass("BM_poisson",
    contains = "scalar_model",
    methods = list(
        initialize = function(membership_type,adj,...)
        {
            .self$initFields(membership_name = membership_type,
                            model_name = "poisson",
                            adj = adj,
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
            if(membership_name=='LBM')
            {
                return(
                    memberships[[Q]]$Z1
                     %*%
                    model_parameters[[Q]]$lambda
                     %*%
                    t(memberships[[Q]]$Z2)
                )
            }
            else
            {
                return(
                    memberships[[Q]]$Z
                     %*%
                    model_parameters[[Q]]$lambda
                     %*%
                    t(memberships[[Q]]$Z)
                )
            }
        }
    )
)

